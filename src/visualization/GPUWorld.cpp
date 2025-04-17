#include "pch.h"
#include "Window.h"
#include "DirectXHelpers.h"
#include "ShaderHelper.h"
#include "GPUWorld.h"
#include "GPUStats.h"

// Constructor
GPUWorld::GPUWorld(std::shared_ptr<Window> pWindow)
    : m_pWindow(pWindow)
    , m_pCamera(std::make_shared<GPUCamera>())
{
    // Set default camera settings
    m_pCamera->setPosition(float3(0.0f, 0.0f, -5.0f));
    m_pCamera->setDirection(float3(0.0f, 0.0f, 1.0f));
    m_pCamera->setFOV(45.0f);
    
    // Get window dimensions
    RECT clientRect;
    GetClientRect(m_pWindow->getWindowHandle(), &clientRect);
    float aspectRatio = static_cast<float>(clientRect.right - clientRect.left) / 
                       static_cast<float>(clientRect.bottom - clientRect.top);
    
    m_pCamera->setAspectRatio(aspectRatio);
    
    // Initialize Direct3D rendering resources
    initializeRenderResources();
}

// Destructor
GPUWorld::~GPUWorld()
{
    // Wait for the GPU to finish all operations before exiting
    Microsoft::WRL::ComPtr<ID3D12Fence> fence;
    ThrowIfFailed(m_pWindow->getDevice()->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&fence)));
    
    auto gpuQueue = m_pWindow->createOrGetGPUQueue();
    uint64_t fenceValue = 1;
    ThrowIfFailed(gpuQueue->getQueue()->Signal(fence.Get(), fenceValue));
    
    // Wait for the fence
    if (fence->GetCompletedValue() < fenceValue)
    {
        HANDLE eventHandle = CreateEvent(nullptr, FALSE, FALSE, nullptr);
        if (eventHandle)
        {
            ThrowIfFailed(fence->SetEventOnCompletion(fenceValue, eventHandle));
            WaitForSingleObject(eventHandle, INFINITE);
            CloseHandle(eventHandle);
        }
    }
}

// Create a new mesh
std::shared_ptr<GPUMesh> GPUWorld::createMesh()
{
    auto device = m_pWindow->getDevice();
    return std::make_shared<GPUMesh>(device);
}

// Add a mesh to the scene
void GPUWorld::addMesh(std::shared_ptr<GPUMesh> mesh)
{
    m_pMeshes.push_back(mesh);
}

// Remove a mesh from the scene
void GPUWorld::removeMesh(std::shared_ptr<GPUMesh> mesh)
{
    auto it = std::find(m_pMeshes.begin(), m_pMeshes.end(), mesh);
    if (it != m_pMeshes.end())
    {
        m_pMeshes.erase(it);
    }
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

// Initialize rendering resources
void GPUWorld::initializeRenderResources()
{
    auto device = m_pWindow->getDevice();
    
    // Create root signature
    CD3DX12_DESCRIPTOR_RANGE ranges[1];
    ranges[0].Init(D3D12_DESCRIPTOR_RANGE_TYPE_CBV, 1, 0); // b0 is Transform buffer
    
    CD3DX12_ROOT_PARAMETER rootParameters[1];
    rootParameters[0].InitAsDescriptorTable(1, &ranges[0], D3D12_SHADER_VISIBILITY_VERTEX);
    
    CD3DX12_ROOT_SIGNATURE_DESC rootSignatureDesc;
    rootSignatureDesc.Init(_countof(rootParameters), rootParameters, 0, nullptr, D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT);
    
    m_rootSignature = CreateRootSignature(device.Get(), rootSignatureDesc);
    
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
        std::wstring shaderPath = L"visualization/Shaders/";
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
    psoDesc.pRootSignature = m_rootSignature.Get();
    
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
    ThrowIfFailed(device->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&m_pipelineState)));
    
    // Create the constant buffer for transformation matrices
    const UINT constantBufferSize = (sizeof(TransformBuffer) + 255) & ~255;
    
    // Create the constant buffer resource
    m_transformBufferResource = CreateBuffer(
        device.Get(),
        static_cast<uint32_t>(constantBufferSize),
        D3D12_RESOURCE_FLAG_NONE,
        D3D12_RESOURCE_STATE_GENERIC_READ);
    
    // Create constant buffer descriptor heap
    D3D12_DESCRIPTOR_HEAP_DESC cbvHeapDesc = {};
    cbvHeapDesc.NumDescriptors = 1;
    cbvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
    cbvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
    ThrowIfFailed(device->CreateDescriptorHeap(&cbvHeapDesc, IID_PPV_ARGS(&m_cbvHeap)));
    
    // Create the constant buffer view
    D3D12_CONSTANT_BUFFER_VIEW_DESC cbvDesc = {};
    cbvDesc.BufferLocation = m_transformBufferResource->GetGPUVirtualAddress();
    cbvDesc.SizeInBytes = constantBufferSize;
    device->CreateConstantBufferView(&cbvDesc, m_cbvHeap->GetCPUDescriptorHandleForHeapStart());
    
    // Map the constant buffer
    CD3DX12_RANGE readRange(0, 0); // We do not intend to read from this resource on the CPU
    ThrowIfFailed(m_transformBufferResource->Map(0, &readRange, reinterpret_cast<void**>(&m_transformBufferData)));
}

// Draw all meshes
void GPUWorld::drawMeshesIntoWindow(GPUStats* pStats)
{
    // Get GPU queue
    auto gpuQueue = m_pWindow->createOrGetGPUQueue();
    
    // Get swap chain
    auto swapChain = m_pWindow->getSwapChain();
    
    // Get current back buffer index
    UINT backBufferIndex = swapChain->GetCurrentBackBufferIndex();
    
    // Get device
    auto device = m_pWindow->getDevice();
    
    // Create render target view heap
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> rtvHeap;
    D3D12_DESCRIPTOR_HEAP_DESC rtvHeapDesc = {};
    rtvHeapDesc.NumDescriptors = 2; // Double buffering
    rtvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;
    rtvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
    ThrowIfFailed(device->CreateDescriptorHeap(&rtvHeapDesc, IID_PPV_ARGS(&rtvHeap)));
    
    // Get rtv descriptor size
    UINT rtvDescriptorSize = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);
    
    // Get back buffers
    Microsoft::WRL::ComPtr<ID3D12Resource> renderTargets[2];
    for (UINT i = 0; i < 2; i++)
    {
        ThrowIfFailed(swapChain->GetBuffer(i, IID_PPV_ARGS(&renderTargets[i])));
        
        CD3DX12_CPU_DESCRIPTOR_HANDLE rtvHandle(
            rtvHeap->GetCPUDescriptorHandleForHeapStart(),
            i,
            rtvDescriptorSize);
        
        device->CreateRenderTargetView(renderTargets[i].Get(), nullptr, rtvHandle);
    }
    
    // Create depth/stencil buffer and view
    Microsoft::WRL::ComPtr<ID3D12Resource> depthStencilBuffer;
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> dsvHeap;
    
    D3D12_DESCRIPTOR_HEAP_DESC dsvHeapDesc = {};
    dsvHeapDesc.NumDescriptors = 1;
    dsvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
    dsvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
    ThrowIfFailed(device->CreateDescriptorHeap(&dsvHeapDesc, IID_PPV_ARGS(&dsvHeap)));
    
    // Get window dimensions
    RECT clientRect;
    GetClientRect(m_pWindow->getWindowHandle(), &clientRect);
    UINT width = clientRect.right - clientRect.left;
    UINT height = clientRect.bottom - clientRect.top;
    
    // Create depth stencil buffer
    D3D12_RESOURCE_DESC depthStencilDesc = {};
    depthStencilDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
    depthStencilDesc.Alignment = 0;
    depthStencilDesc.Width = width;
    depthStencilDesc.Height = height;
    depthStencilDesc.DepthOrArraySize = 1;
    depthStencilDesc.MipLevels = 1;
    depthStencilDesc.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
    depthStencilDesc.SampleDesc.Count = 1;
    depthStencilDesc.SampleDesc.Quality = 0;
    depthStencilDesc.Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN;
    depthStencilDesc.Flags = D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL;
    
    D3D12_CLEAR_VALUE clearValue = {};
    clearValue.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
    clearValue.DepthStencil.Depth = 1.0f;
    clearValue.DepthStencil.Stencil = 0;
    
    CD3DX12_HEAP_PROPERTIES heapProperties(D3D12_HEAP_TYPE_DEFAULT);
    
    ThrowIfFailed(device->CreateCommittedResource(
        &heapProperties,
        D3D12_HEAP_FLAG_NONE,
        &depthStencilDesc,
        D3D12_RESOURCE_STATE_DEPTH_WRITE,
        &clearValue,
        IID_PPV_ARGS(&depthStencilBuffer)));
    
    // Create depth stencil view
    D3D12_DEPTH_STENCIL_VIEW_DESC dsvDesc = {};
    dsvDesc.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
    dsvDesc.ViewDimension = D3D12_DSV_DIMENSION_TEXTURE2D;
    dsvDesc.Flags = D3D12_DSV_FLAG_NONE;
    
    device->CreateDepthStencilView(
        depthStencilBuffer.Get(),
        &dsvDesc,
        dsvHeap->GetCPUDescriptorHandleForHeapStart());
    
    // Begin command list recording
    auto commandList = gpuQueue->beginRecording();
    if (pStats)
    {
        pStats->begin(*commandList.Get());
    }
    
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
    
    commandList->RSSetViewports(1, &viewport);
    commandList->RSSetScissorRects(1, &scissorRect);
    
    // Transition render target to render target state
    CD3DX12_RESOURCE_BARRIER renderTargetBarrier = 
        CD3DX12_RESOURCE_BARRIER::Transition(
            renderTargets[backBufferIndex].Get(),
            D3D12_RESOURCE_STATE_PRESENT,
            D3D12_RESOURCE_STATE_RENDER_TARGET);
    
    commandList->ResourceBarrier(1, &renderTargetBarrier);
    
    // Get render target view handle
    CD3DX12_CPU_DESCRIPTOR_HANDLE rtvHandle(
        rtvHeap->GetCPUDescriptorHandleForHeapStart(),
        backBufferIndex,
        rtvDescriptorSize);
    
    // Get depth stencil view handle
    CD3DX12_CPU_DESCRIPTOR_HANDLE dsvHandle(dsvHeap->GetCPUDescriptorHandleForHeapStart());
    
    // Set render targets
    commandList->OMSetRenderTargets(1, &rtvHandle, FALSE, &dsvHandle);
    
    // Clear render target and depth stencil
    const float clearColor[] = { 0.0f, 0.2f, 0.4f, 1.0f }; // Dark blue background
    commandList->ClearRenderTargetView(rtvHandle, clearColor, 0, nullptr);
    commandList->ClearDepthStencilView(dsvHandle, D3D12_CLEAR_FLAG_DEPTH | D3D12_CLEAR_FLAG_STENCIL, 1.0f, 0, 0, nullptr);
    
    // Update the transform constant buffer
    if (m_transformBufferData)
    {
        // Update matrices
        m_transformBufferMatrix.World = DirectX::XMMatrixIdentity();
        m_transformBufferMatrix.View = m_pCamera->getViewMatrix();
        m_transformBufferMatrix.Projection = m_pCamera->getProjectionMatrix();
        
        // Copy to GPU memory
        memcpy(m_transformBufferData, &m_transformBufferMatrix, sizeof(TransformBuffer));
    }
    
    // Set the descriptor heap, root signature, and pipeline state
    ID3D12DescriptorHeap* heaps[] = { m_cbvHeap.Get() };
    commandList->SetDescriptorHeaps(_countof(heaps), heaps);
    commandList->SetGraphicsRootSignature(m_rootSignature.Get());
    commandList->SetPipelineState(m_pipelineState.Get());
    commandList->SetGraphicsRootDescriptorTable(0, m_cbvHeap->GetGPUDescriptorHandleForHeapStart());
    
    // Draw meshes
    for (const auto& mesh : m_pMeshes)
    {
        // Set primitive topology
        commandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        
        // Set vertex buffer
        D3D12_VERTEX_BUFFER_VIEW vbv = mesh->getVertexBufferView();
        commandList->IASetVertexBuffers(0, 1, &vbv);
        
        // Set index buffer
        D3D12_INDEX_BUFFER_VIEW ibv = mesh->getIndexBufferView();
        commandList->IASetIndexBuffer(&ibv);
        
        // Draw
        commandList->DrawIndexedInstanced(mesh->getIndexCount(), 1, 0, 0, 0);
    }
    
    // Transition render target back to present state
    renderTargetBarrier = CD3DX12_RESOURCE_BARRIER::Transition(
        renderTargets[backBufferIndex].Get(),
        D3D12_RESOURCE_STATE_RENDER_TARGET,
        D3D12_RESOURCE_STATE_PRESENT);
    
    commandList->ResourceBarrier(1, &renderTargetBarrier);
    
    // Execute command list
    if (pStats)
    {
        pStats->end(*commandList.Get());
    }
    gpuQueue->execute(commandList);
    
    // Present the frame
    ThrowIfFailed(swapChain->Present(1, 0));
} 