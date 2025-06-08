#include "pch.h"
#include "GPUText.h"
#include "GPUQueue.h"
#include "SwapChain.h"
#include "DirectXHelpers.h"
#include <vector>
#include <sstream>

GPUText::GPUText(std::shared_ptr<GPUFont> pFont)
    : m_pFont(pFont)
    , m_vLeftTop(0.0f, 0.0f)
{
}

void GPUText::setLeftTop(const float2& vLeftTop)
{
    m_vLeftTop = vLeftTop;
}

int GPUText::printf(const char* format, ...)
{
    if (!format) {
        return 0;
    }
    
    // Handle variable arguments
    va_list args;
    va_start(args, format);
    
    // First, determine the required buffer size
    va_list args_copy;
    va_copy(args_copy, args);
    int required_size = vsnprintf(nullptr, 0, format, args_copy);
    va_end(args_copy);
    
    if (required_size < 0) {
        va_end(args);
        return 0;
    }
    
    // Allocate buffer and format the string
    std::vector<char> buffer(required_size + 1);
    int formatted_chars = vsnprintf(buffer.data(), buffer.size(), format, args);
    va_end(args);
    
    if (formatted_chars < 0) {
        return 0;
    }
    
    // Convert to string and split by newlines
    std::string formatted_text(buffer.data());
    
    // Clear existing lines
    m_sLines.clear();
    
    // Split the formatted text by newlines
    std::istringstream stream(formatted_text);
    std::string line;
    
    while (std::getline(stream, line)) {
        m_sLines.push_back(line);
    }
    
    // If the formatted text doesn't end with a newline and we have content,
    // make sure we don't lose the last line
    if (!formatted_text.empty() && formatted_text.back() != '\n' && m_sLines.empty()) {
        m_sLines.push_back(formatted_text);
    }
    
    return formatted_chars;
}

void GPUText::generateTextQuads(std::vector<TextVertex>& vertices, std::vector<uint16_t>& indices, 
                               float2 screenSize)
{
    vertices.clear();
    indices.clear();
    
    if (m_sLines.empty() || !m_pFont) {
        return;
    }
    
    float currentX = m_vLeftTop.x;
    float currentY = m_vLeftTop.y;
    float lineHeight = m_pFont->getLineHeight();
    
    uint16_t vertexIndex = 0;
    
    for (const auto& line : m_sLines) {
        currentX = m_vLeftTop.x; // Reset X for each line
        
        for (char c : line) {
            const GlyphInfo* glyphInfo = m_pFont->getGlyphInfo(c);
            if (!glyphInfo) {
                continue; // Skip unknown characters
            }
            
            // Calculate glyph position
            float glyphX = currentX + glyphInfo->bearing.x;
            float glyphY = currentY + lineHeight + glyphInfo->bearing.y;
            float glyphWidth = glyphInfo->size.x;
            float glyphHeight = glyphInfo->size.y;

            // Create quad vertices (2 triangles)
            TextVertex quad[4];
            
            // Top-left
            quad[0].position = float2(glyphX, glyphY);
            quad[0].texCoord = glyphInfo->texCoords[0];
            
            // Top-right
            quad[1].position = float2(glyphX + glyphWidth, glyphY);
            quad[1].texCoord = float2(glyphInfo->texCoords[1].x, glyphInfo->texCoords[0].y);
            
            // Bottom-right
            quad[2].position = float2(glyphX + glyphWidth, glyphY + glyphHeight);
            quad[2].texCoord = glyphInfo->texCoords[1];
            
            // Bottom-left
            quad[3].position = float2(glyphX, glyphY + glyphHeight);
            quad[3].texCoord = float2(glyphInfo->texCoords[0].x, glyphInfo->texCoords[1].y);
            
            // Add vertices
            for (int i = 0; i < 4; i++) {
                vertices.push_back(quad[i]);
            }
            
            // Add indices for two triangles (CCW winding)
            indices.push_back(vertexIndex + 0); // Triangle 1
            indices.push_back(vertexIndex + 1);
            indices.push_back(vertexIndex + 2);
            
            indices.push_back(vertexIndex + 0); // Triangle 2
            indices.push_back(vertexIndex + 2);
            indices.push_back(vertexIndex + 3);
            
            vertexIndex += 4;
            
            // Advance to next character position
            currentX += glyphInfo->advance;
        }
        
        // Move to next line
        currentY += lineHeight;
    }
}

void GPUText::updateVertexBuffer(const std::vector<TextVertex>& vertices, 
                                const std::vector<uint16_t>& indices, 
                                ID3D12Device* pDevice)
{
    if (vertices.empty() || indices.empty()) {
        return;
    }
    
    // Calculate buffer sizes
    UINT vertexBufferSize = static_cast<UINT>(vertices.size() * sizeof(TextVertex));
    UINT indexBufferSize = static_cast<UINT>(indices.size() * sizeof(uint16_t));
    
    // Create or recreate vertex buffer if size changed
    if (!m_vertexBuffer || m_vertexCount != vertices.size()) {
        m_vertexBuffer = CreateBuffer(pDevice, vertexBufferSize, 
                                     D3D12_RESOURCE_FLAG_NONE, 
                                     D3D12_RESOURCE_STATE_GENERIC_READ);
        
        m_vertexBufferView.BufferLocation = m_vertexBuffer->GetGPUVirtualAddress();
        m_vertexBufferView.StrideInBytes = sizeof(TextVertex);
        m_vertexBufferView.SizeInBytes = vertexBufferSize;
        m_vertexCount = static_cast<uint32_t>(vertices.size());
    }
    
    // Create or recreate index buffer if size changed
    if (!m_indexBuffer || m_indexCount != indices.size()) {
        m_indexBuffer = CreateBuffer(pDevice, indexBufferSize,
                                    D3D12_RESOURCE_FLAG_NONE,
                                    D3D12_RESOURCE_STATE_GENERIC_READ);
        
        m_indexBufferView.BufferLocation = m_indexBuffer->GetGPUVirtualAddress();
        m_indexBufferView.Format = DXGI_FORMAT_R16_UINT;
        m_indexBufferView.SizeInBytes = indexBufferSize;
        m_indexCount = static_cast<uint32_t>(indices.size());
    }
    
    // Upload vertex data
    UploadToBuffer(m_vertexBuffer.Get(), vertices.data(), static_cast<uint32_t>(vertices.size()));
    
    // Upload index data
    UploadToBuffer(m_indexBuffer.Get(), indices.data(), static_cast<uint32_t>(indices.size()));
}

void GPUText::updateConstantBuffer(const float2& screenSize, ID3D12Device* pDevice)
{
    
    // Create constant buffer if it doesn't exist
    if (!m_constantBuffer) {
        const UINT constantBufferSize = (sizeof(TextParams) + 255) & ~255; // 256-byte aligned
        
        m_constantBuffer = CreateBuffer(pDevice, constantBufferSize,
                                       D3D12_RESOURCE_FLAG_NONE,
                                       D3D12_RESOURCE_STATE_GENERIC_READ);
        
        // Map the constant buffer for persistent writing
        CD3DX12_RANGE readRange(0, 0);
        ThrowIfFailed(m_constantBuffer->Map(0, &readRange, reinterpret_cast<void**>(&m_constantBufferData)));
    }
    
    // Update constant buffer data
    if (m_constantBufferData) {
        TextParams params;
        params.textColor = m_textColor;
        params.screenSize = screenSize;
        params.padding = float2(0.0f, 0.0f);
        
        memcpy(m_constantBufferData, &params, sizeof(TextParams));
    }
}

void GPUText::ensureDescriptorHeaps(ID3D12Device* pDevice)
{
    // Create combined CBV/SRV descriptor heap
    if (!m_descriptorHeap) {
        D3D12_DESCRIPTOR_HEAP_DESC heapDesc = {};
        heapDesc.NumDescriptors = 2; // 0: CBV for text params, 1: SRV for font texture
        heapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
        heapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
        ThrowIfFailed(pDevice->CreateDescriptorHeap(&heapDesc, IID_PPV_ARGS(&m_descriptorHeap)));
        
        UINT descriptorSize = pDevice->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);
        CD3DX12_CPU_DESCRIPTOR_HANDLE heapStart(m_descriptorHeap->GetCPUDescriptorHandleForHeapStart());
        
        // Create CBV for text parameters at index 0
        if (m_constantBuffer) {
            D3D12_CONSTANT_BUFFER_VIEW_DESC cbvDesc = {};
            cbvDesc.BufferLocation = m_constantBuffer->GetGPUVirtualAddress();
            cbvDesc.SizeInBytes = (sizeof(TextParams) + 255) & ~255;
            pDevice->CreateConstantBufferView(&cbvDesc, heapStart);
        }
        
        // Create SRV for font atlas at index 1
        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
        srvDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
        srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
        srvDesc.Texture2D.MipLevels = 1;
        
        auto fontResource = m_pFont->getResource();
        if (fontResource && fontResource->getResource()) {
            CD3DX12_CPU_DESCRIPTOR_HANDLE srvHandle(heapStart, 1, descriptorSize);
            pDevice->CreateShaderResourceView(fontResource->getResource().Get(), &srvDesc, srvHandle);
        }
    }
}

void GPUText::render(SwapChain* pSwapChain, ID3D12RootSignature* pSharedRS,
                    ID3D12GraphicsCommandList* pCmdList)
{
    if (m_sLines.empty() || !m_pFont || !pSharedRS || !pCmdList || !pSwapChain) {
        return;
    }
    
    // Get swap chain dimensions
    DXGI_SWAP_CHAIN_DESC1 swapChainDesc;
    ThrowIfFailed(pSwapChain->getSwapChain()->GetDesc1(&swapChainDesc));
    float2 screenSize = float2(static_cast<float>(swapChainDesc.Width), static_cast<float>(swapChainDesc.Height));
    
    // Get device from command list
    Microsoft::WRL::ComPtr<ID3D12Device> pDevice;
    pCmdList->GetDevice(IID_PPV_ARGS(&pDevice));
    
    // Generate text geometry
    std::vector<TextVertex> vertices;
    std::vector<uint16_t> indices;
    generateTextQuads(vertices, indices, screenSize);
    
    if (vertices.empty() || indices.empty()) {
        return; // Nothing to render
    }
    
    // Update buffers
    updateVertexBuffer(vertices, indices, pDevice.Get());
    updateConstantBuffer(screenSize, pDevice.Get());
    
    // Ensure descriptor heaps exist
    ensureDescriptorHeaps(pDevice.Get());
    
    // Update CBV if constant buffer was just created
    if (m_constantBuffer && m_descriptorHeap) {
        UINT descriptorSize = pDevice->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);
        CD3DX12_CPU_DESCRIPTOR_HANDLE heapStart(m_descriptorHeap->GetCPUDescriptorHandleForHeapStart());
        
        D3D12_CONSTANT_BUFFER_VIEW_DESC cbvDesc = {};
        cbvDesc.BufferLocation = m_constantBuffer->GetGPUVirtualAddress();
        cbvDesc.SizeInBytes = (sizeof(TextParams) + 255) & ~255;
        pDevice->CreateConstantBufferView(&cbvDesc, heapStart);
    }
    
    // Get back buffer resource and RTV handle from SwapChain
    ID3D12Resource* backBufferResource = pSwapChain->getBBColor();
    D3D12_CPU_DESCRIPTOR_HANDLE rtvHandle = pSwapChain->getBBColorCPUHandle();
    D3D12_CPU_DESCRIPTOR_HANDLE dsvHandle = pSwapChain->getBBDepthCPUHandle();
    
    // Transition back buffer to render target state
    CD3DX12_RESOURCE_BARRIER renderTargetBarrier = 
        CD3DX12_RESOURCE_BARRIER::Transition(
            backBufferResource,
            D3D12_RESOURCE_STATE_PRESENT,
            D3D12_RESOURCE_STATE_RENDER_TARGET);
    pCmdList->ResourceBarrier(1, &renderTargetBarrier);
    
    // Set render targets with depth stencil view
    pCmdList->OMSetRenderTargets(1, &rtvHandle, FALSE, &dsvHandle);
    
    // Set viewport and scissor rect
    D3D12_VIEWPORT viewport = {};
    viewport.TopLeftX = 0.0f;
    viewport.TopLeftY = 0.0f;
    viewport.Width = screenSize.x;
    viewport.Height = screenSize.y;
    viewport.MinDepth = 0.0f;
    viewport.MaxDepth = 1.0f;
    pCmdList->RSSetViewports(1, &viewport);
    
    D3D12_RECT scissorRect = {};
    scissorRect.left = 0;
    scissorRect.top = 0;
    scissorRect.right = static_cast<LONG>(screenSize.x);
    scissorRect.bottom = static_cast<LONG>(screenSize.y);
    pCmdList->RSSetScissorRects(1, &scissorRect);
    
    // Get text PSO from font
    auto textPSO = m_pFont->getTextPSO(pSharedRS);
    if (!textPSO) {
        return;
    }
    
    // Set pipeline state and root signature
    pCmdList->SetGraphicsRootSignature(pSharedRS);
    pCmdList->SetPipelineState(textPSO);
    
    // Set primitive topology
    pCmdList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    
    // Set vertex and index buffers
    pCmdList->IASetVertexBuffers(0, 1, &m_vertexBufferView);
    pCmdList->IASetIndexBuffer(&m_indexBufferView);
    
    // Set descriptor heaps
    ID3D12DescriptorHeap* heaps[] = { m_descriptorHeap.Get() };
    pCmdList->SetDescriptorHeaps(_countof(heaps), heaps);
    
    // Set root signature parameters
    UINT descriptorSize = pDevice->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);
    CD3DX12_GPU_DESCRIPTOR_HANDLE heapStart(m_descriptorHeap->GetGPUDescriptorHandleForHeapStart());
    
    // Bind constant buffer (b0 for text parameters)
    pCmdList->SetGraphicsRootDescriptorTable(0, heapStart); // CBV at index 0
    
    // Bind font texture (t0)
    pCmdList->SetGraphicsRootDescriptorTable(1, CD3DX12_GPU_DESCRIPTOR_HANDLE(heapStart, 1, descriptorSize)); // SRV at index 1
    
    // Draw text
    pCmdList->DrawIndexedInstanced(m_indexCount, 1, 0, 0, 0);
    
    // Transition back buffer back to present state
    renderTargetBarrier = CD3DX12_RESOURCE_BARRIER::Transition(
        backBufferResource,
        D3D12_RESOURCE_STATE_RENDER_TARGET,
        D3D12_RESOURCE_STATE_PRESENT);
    pCmdList->ResourceBarrier(1, &renderTargetBarrier);
}
