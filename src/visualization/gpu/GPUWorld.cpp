#include "pch.h"
#include "Window.h"
#include "DirectXHelpers.h"
#include "ShaderHelper.h"
#include "GPUWorld.h"
#include "GPUStats.h"
#include "GPUMesh.h"
#include "IObjectVis.h"

// Constructor
GPUWorld::GPUWorld(std::shared_ptr<Window> pWindow, GPUQueue* pGpuQueue)
    : m_pWindow(pWindow)
    , m_pCamera(std::make_shared<GPUCamera>())
{
    // Set default camera settings
    m_pCamera->setPosition(float3(0.0f, 0.0f, -5.0f));
    m_pCamera->setDirection(float3(0.0f, 0.0f, 1.0f));
    m_pCamera->setFOV(45.0f);

    // Create GPU font using the temporary queue
    m_pFont = std::make_shared<GPUFont>(48, pGpuQueue); // 48px font size
    
    // Initialize Direct3D rendering resources
    initializeRenderResources();
}

// Destructor
GPUWorld::~GPUWorld()
{
    // Resources will be automatically released by ComPtr destructors
    // GPU synchronization should be handled by the caller when needed
}

// Add an object to the scene
void GPUWorld::addObject(std::weak_ptr<IObjectVis> pObject)
{
    m_pObjects.push_back(pObject);
}

// Get current camera
std::shared_ptr<GPUCamera> GPUWorld::getCamera()
{
    return m_pCamera;
}

// Set a new camera
void GPUWorld::setCamera(std::shared_ptr<GPUCamera> camera)
{
    if (camera)
    {
        m_pCamera = camera;
    }
}

// Get current font
std::shared_ptr<GPUFont> GPUWorld::getFont()
{
    return m_pFont;
}

// Initialize rendering resources
void GPUWorld::initializeRenderResources()
{
    auto pDevice = m_pWindow->getDevice();
    
    // Create shared root signature for both mesh and text rendering
    CD3DX12_DESCRIPTOR_RANGE ranges[2];
    ranges[0].Init(D3D12_DESCRIPTOR_RANGE_TYPE_CBV, 2, 0);  // b0-b1: Constant buffers (view/projection + text params)
    ranges[1].Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, 4, 0);  // t0-t3: Textures (font atlas, future textures)
    
    CD3DX12_ROOT_PARAMETER rootParameters[3];
    rootParameters[0].InitAsDescriptorTable(1, &ranges[0], D3D12_SHADER_VISIBILITY_ALL);
    rootParameters[1].InitAsDescriptorTable(1, &ranges[1], D3D12_SHADER_VISIBILITY_ALL);
    rootParameters[2].InitAsConstants(16, 2, 0, D3D12_SHADER_VISIBILITY_VERTEX);  // 16 floats for 4x4 world matrix at b2
    
    // Static sampler for texture sampling (font atlas, etc.)
    D3D12_STATIC_SAMPLER_DESC sampler = {};
    sampler.Filter = D3D12_FILTER_MIN_MAG_LINEAR_MIP_POINT;
    sampler.AddressU = D3D12_TEXTURE_ADDRESS_MODE_CLAMP;
    sampler.AddressV = D3D12_TEXTURE_ADDRESS_MODE_CLAMP;
    sampler.AddressW = D3D12_TEXTURE_ADDRESS_MODE_CLAMP;
    sampler.MipLODBias = 0;
    sampler.MaxAnisotropy = 0;
    sampler.ComparisonFunc = D3D12_COMPARISON_FUNC_NEVER;
    sampler.BorderColor = D3D12_STATIC_BORDER_COLOR_TRANSPARENT_BLACK;
    sampler.MinLOD = 0.0f;
    sampler.MaxLOD = D3D12_FLOAT32_MAX;
    sampler.ShaderRegister = 0;
    sampler.RegisterSpace = 0;
    sampler.ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL;
    
    CD3DX12_ROOT_SIGNATURE_DESC rootSignatureDesc;
    rootSignatureDesc.Init(_countof(rootParameters), rootParameters, 1, &sampler, D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT);
    
    m_pRootSignature = CreateRootSignature(pDevice, rootSignatureDesc);
    
    // Create the pipeline state, which includes compiling and loading shaders
    
    // Load and compile shaders
    Microsoft::WRL::ComPtr<ID3DBlob> vertexShader;
    Microsoft::WRL::ComPtr<ID3DBlob> pixelShader;
    
#if defined(_DEBUG)
    // Enable better shader debugging with the graphics debugging tools
    UINT compileFlags = D3DCOMPILE_DEBUG | D3DCOMPILE_SKIP_OPTIMIZATION;
#else
    UINT compileFlags = 0;
#endif
    
    // Use our shader helper to load/compile shaders
    ShaderHelper& shaderHelper = ShaderHelper::getInstance();
    
    {
        // Try to load pre-compiled shaders
        std::wstring shaderPath = L"Shaders/";
        vertexShader = shaderHelper.loadCompiledShader(shaderPath + L"VertexShader.cso");
        pixelShader = shaderHelper.loadCompiledShader(shaderPath + L"PixelShader.cso");
    }

    {
        // Fallback: Compile shaders at runtime
        std::wstring shaderPath = L"visualization/gpu/Shaders/";
        if (!vertexShader)
            vertexShader = shaderHelper.loadShader(shaderPath + L"VertexShader.hlsl", "main", "vs_5_0", compileFlags);
        if (!pixelShader)
            pixelShader = shaderHelper.loadShader(shaderPath + L"PixelShader.hlsl", "main", "ps_5_0", compileFlags);
    }

    // Define the vertex input layout
    D3D12_INPUT_ELEMENT_DESC inputElementDescs[] =
    {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 }
    };
    
    // Setup the pipeline state using the traditional, non-stream API for compatibility
    D3D12_GRAPHICS_PIPELINE_STATE_DESC psoDesc = {};
    
    // Setup the root signature
    psoDesc.pRootSignature = m_pRootSignature.Get();
    
    // Setup the vertex shader
    psoDesc.VS.pShaderBytecode = vertexShader->GetBufferPointer();
    psoDesc.VS.BytecodeLength = vertexShader->GetBufferSize();
    
    // Setup the pixel shader
    psoDesc.PS.pShaderBytecode = pixelShader->GetBufferPointer();
    psoDesc.PS.BytecodeLength = pixelShader->GetBufferSize();
    
    // Setup rasterizer state
    psoDesc.RasterizerState.FillMode = D3D12_FILL_MODE_WIREFRAME;
    psoDesc.RasterizerState.CullMode = D3D12_CULL_MODE_BACK;
    psoDesc.RasterizerState.FrontCounterClockwise = FALSE;
    psoDesc.RasterizerState.DepthBias = D3D12_DEFAULT_DEPTH_BIAS;
    psoDesc.RasterizerState.DepthBiasClamp = D3D12_DEFAULT_DEPTH_BIAS_CLAMP;
    psoDesc.RasterizerState.SlopeScaledDepthBias = D3D12_DEFAULT_SLOPE_SCALED_DEPTH_BIAS;
    psoDesc.RasterizerState.DepthClipEnable = FALSE; // TODO: figure out why I see nothing with 'TRUE' here
    psoDesc.RasterizerState.MultisampleEnable = FALSE;
    psoDesc.RasterizerState.AntialiasedLineEnable = FALSE;
    psoDesc.RasterizerState.ForcedSampleCount = 0;
    psoDesc.RasterizerState.ConservativeRaster = D3D12_CONSERVATIVE_RASTERIZATION_MODE_OFF;
    
    // Setup blend state
    D3D12_BLEND_DESC blendDesc = {};
    blendDesc.AlphaToCoverageEnable = FALSE;
    blendDesc.IndependentBlendEnable = FALSE;
    blendDesc.RenderTarget[0].BlendEnable = FALSE;
    blendDesc.RenderTarget[0].LogicOpEnable = FALSE;
    blendDesc.RenderTarget[0].SrcBlend = D3D12_BLEND_ONE;
    blendDesc.RenderTarget[0].DestBlend = D3D12_BLEND_ZERO;
    blendDesc.RenderTarget[0].BlendOp = D3D12_BLEND_OP_ADD;
    blendDesc.RenderTarget[0].SrcBlendAlpha = D3D12_BLEND_ONE;
    blendDesc.RenderTarget[0].DestBlendAlpha = D3D12_BLEND_ZERO;
    blendDesc.RenderTarget[0].BlendOpAlpha = D3D12_BLEND_OP_ADD;
    blendDesc.RenderTarget[0].LogicOp = D3D12_LOGIC_OP_NOOP;
    blendDesc.RenderTarget[0].RenderTargetWriteMask = D3D12_COLOR_WRITE_ENABLE_ALL;
    psoDesc.BlendState = blendDesc;
    
    // Setup depth stencil state
    D3D12_DEPTH_STENCIL_DESC depthStencilDesc = {};
    depthStencilDesc.DepthEnable = TRUE;
    depthStencilDesc.DepthWriteMask = D3D12_DEPTH_WRITE_MASK_ALL;
    depthStencilDesc.DepthFunc = D3D12_COMPARISON_FUNC_LESS;
    depthStencilDesc.StencilEnable = FALSE;
    psoDesc.DepthStencilState = depthStencilDesc;
    
    // Setup input layout
    psoDesc.InputLayout = { inputElementDescs, _countof(inputElementDescs) };
    
    // Setup primitive topology
    psoDesc.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
    
    // Setup render targets
    psoDesc.NumRenderTargets = 1;
    psoDesc.RTVFormats[0] = DXGI_FORMAT_R8G8B8A8_UNORM;
    psoDesc.DSVFormat = DXGI_FORMAT_D24_UNORM_S8_UINT;
    
    // Setup MSAA
    psoDesc.SampleDesc.Count = 1;
    psoDesc.SampleDesc.Quality = 0;
    psoDesc.SampleMask = UINT_MAX;  // Enable all samples
    
    // Setup misc flags
    psoDesc.NodeMask = 0;
    psoDesc.CachedPSO.pCachedBlob = nullptr;
    psoDesc.CachedPSO.CachedBlobSizeInBytes = 0;
    psoDesc.Flags = D3D12_PIPELINE_STATE_FLAG_NONE;
    
    // Create the pipeline state
    ThrowIfFailed(pDevice->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&m_pPipelineState)));
    
    // Create the constant buffer for transformation matrices
    const UINT constantBufferSize = (sizeof(TransformBuffer) + 255) & ~255;
    
    // Create the constant buffer resource
    m_pTransformRes = CreateBuffer(
        pDevice,
        static_cast<uint32_t>(constantBufferSize),
        D3D12_RESOURCE_FLAG_NONE,
        D3D12_RESOURCE_STATE_GENERIC_READ);
    
    // Create constant buffer descriptor heap
    D3D12_DESCRIPTOR_HEAP_DESC cbvHeapDesc = {};
    cbvHeapDesc.NumDescriptors = 1;
    cbvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
    cbvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
    ThrowIfFailed(pDevice->CreateDescriptorHeap(&cbvHeapDesc, IID_PPV_ARGS(&m_pCBVHeap)));
    
    // Create the constant buffer view
    D3D12_CONSTANT_BUFFER_VIEW_DESC cbvDesc = {};
    cbvDesc.BufferLocation = m_pTransformRes->GetGPUVirtualAddress();
    cbvDesc.SizeInBytes = constantBufferSize;
    pDevice->CreateConstantBufferView(&cbvDesc, m_pCBVHeap->GetCPUDescriptorHandleForHeapStart());
    
    // Map the constant buffer
    CD3DX12_RANGE readRange(0, 0); // We do not intend to read from this resource on the CPU
    ThrowIfFailed(m_pTransformRes->Map(0, &readRange, reinterpret_cast<void**>(&m_pTransformData)));
}

// Draw all meshes and return combined bounding box
box3 GPUWorld::render(SwapChain* pSwapChain, ID3D12GraphicsCommandList* pCmdList)
{
    // Get swap chain dimensions
    DXGI_SWAP_CHAIN_DESC1 swapChainDesc;
    ThrowIfFailed(pSwapChain->getSwapChain()->GetDesc1(&swapChainDesc));
    UINT width = swapChainDesc.Width;
    UINT height = swapChainDesc.Height;

    // Set viewport and scissor rect
    D3D12_VIEWPORT viewport = {};
    viewport.TopLeftX = 0;
    viewport.TopLeftY = 0;
    viewport.Width = static_cast<float>(width);
    viewport.Height = static_cast<float>(height);
    viewport.MinDepth = 0.0f;
    viewport.MaxDepth = 1.0f;
    
    D3D12_RECT scissorRect = {};
    scissorRect.left = 0;
    scissorRect.top = 0;
    scissorRect.right = width;
    scissorRect.bottom = height;
    
    pCmdList->RSSetViewports(1, &viewport);
    pCmdList->RSSetScissorRects(1, &scissorRect);
    
    // Get back buffer resource and handles from SwapChain
    ID3D12Resource* backBufferResource = pSwapChain->getBBColor();
    D3D12_CPU_DESCRIPTOR_HANDLE rtvHandle = pSwapChain->getBBColorCPUHandle();
    D3D12_CPU_DESCRIPTOR_HANDLE dsvHandle = pSwapChain->getBBDepthCPUHandle();
    
    // Transition render target to render target state
    CD3DX12_RESOURCE_BARRIER renderTargetBarrier = 
        CD3DX12_RESOURCE_BARRIER::Transition(
            backBufferResource,
            D3D12_RESOURCE_STATE_PRESENT,
            D3D12_RESOURCE_STATE_RENDER_TARGET);
    
    pCmdList->ResourceBarrier(1, &renderTargetBarrier);
    
    // Set render targets
    pCmdList->OMSetRenderTargets(1, &rtvHandle, FALSE, &dsvHandle);
    
    // Clear render target and depth stencil
    const float clearColor[] = { 0.0f, 0.2f, 0.4f, 1.0f }; // Dark blue background
    pCmdList->ClearRenderTargetView(rtvHandle, clearColor, 0, nullptr);
    pCmdList->ClearDepthStencilView(dsvHandle, D3D12_CLEAR_FLAG_DEPTH | D3D12_CLEAR_FLAG_STENCIL, 1.0f, 0, 0, nullptr);
    
    // Update the transform constant buffer with camera matrices (shared by all meshes)
    if (m_pTransformData)
    {
        m_pCamera->setAspectRatio(width / (float)height);
        // Update view and projection matrices once
        m_transformBufferMatrix.View = m_pCamera->getViewMatrix();
        m_transformBufferMatrix.Projection = m_pCamera->getProjectionMatrix();
        
        // Copy updated transform buffer to GPU memory (View and Projection only)
        memcpy(m_pTransformData, &m_transformBufferMatrix, sizeof(TransformBuffer));
    }
    
    // Set the descriptor heap, root signature, and pipeline state
    ID3D12DescriptorHeap* heaps[] = { m_pCBVHeap.Get() };
    pCmdList->SetDescriptorHeaps(_countof(heaps), heaps);
    pCmdList->SetGraphicsRootSignature(m_pRootSignature.Get());
    pCmdList->SetPipelineState(m_pPipelineState.Get());
    pCmdList->SetGraphicsRootDescriptorTable(0, m_pCBVHeap->GetGPUDescriptorHandleForHeapStart());
    
    // Initialize bounding box accumulator (start with empty box)
    box3 sceneBoundingBox;
    sceneBoundingBox.m_mins = float3(FLT_MAX, FLT_MAX, FLT_MAX);
    sceneBoundingBox.m_maxs = float3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    bool hasValidBounds = false;
    
    // Draw objects - get mesh from each object and render
    for (auto itObject = m_pObjects.begin(); itObject != m_pObjects.end(); )
    {
        auto pObject = itObject->lock();
        if (!pObject)
        {
            itObject = m_pObjects.erase(itObject);
            continue;
        }
        
        // Get the mesh from the object
        auto pMesh = pObject->updateAndGetGpuMesh();
        if (!pMesh)
        {
            ++itObject; // Move to next object if no mesh
            continue; // Skip objects that don't have a mesh
        }
        
        // Set the world matrix as root constants for this specific mesh
        DirectX::XMMATRIX worldMatrix = pMesh->getWorldMatrix();
        pCmdList->SetGraphicsRoot32BitConstants(2, 16, &worldMatrix, 0);  // Root parameter 2, 16 floats (4x4 matrix)
        
        // Set primitive topology
        pCmdList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        
        // Set vertex buffer
        D3D12_VERTEX_BUFFER_VIEW vbv = pMesh->getVertexBufferView();
        pCmdList->IASetVertexBuffers(0, 1, &vbv);
        
        // Set index buffer
        D3D12_INDEX_BUFFER_VIEW ibv = pMesh->getIndexBufferView();
        pCmdList->IASetIndexBuffer(&ibv);
        
        // Draw
        pCmdList->DrawIndexedInstanced(pMesh->getIndexCount(), 1, 0, 0, 0);
        
        // Accumulate bounding box - transform local mesh bounds to world space
        const box3& localBounds = pMesh->getBoundingBox();
        if (!localBounds.isempty())
        {
            // Transform bounding box to world space using mesh transform
            box3 worldBounds = localBounds * pMesh->getTransform();
            
            // Union with scene bounding box
            if (hasValidBounds)
            {
                sceneBoundingBox = sceneBoundingBox | worldBounds;
            }
            else
            {
                sceneBoundingBox = worldBounds;
                hasValidBounds = true;
            }
        }
        
        ++itObject; // Move to next object
    }
    
    // Transition render target back to present state
    renderTargetBarrier = CD3DX12_RESOURCE_BARRIER::Transition(
        backBufferResource,
        D3D12_RESOURCE_STATE_RENDER_TARGET,
        D3D12_RESOURCE_STATE_PRESENT);
    
    pCmdList->ResourceBarrier(1, &renderTargetBarrier);
    
    // Return the computed scene bounding box (empty if no valid meshes)
    if (!hasValidBounds)
    {
        sceneBoundingBox.m_mins = float3(0.0f, 0.0f, 0.0f);
        sceneBoundingBox.m_maxs = float3(0.0f, 0.0f, 0.0f);
    }
    
    return sceneBoundingBox;
} 