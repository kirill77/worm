#pragma once

#include <memory>
#include <cstdarg>
#include <vector>
#include "GPUFont.h"
#include "math/vector.h"
#include <d3d12.h>
#include <wrl/client.h>

// Forward declarations
struct GPUQueue;
struct SwapChain;

struct GPUText
{
public:
    GPUText(std::shared_ptr<GPUFont> pFont);

    void setLeftTop(const float2& vLeftTop);
    int printf(const char* format, ...);
    void render(SwapChain* pSwapChain, ID3D12RootSignature *pSharedRS,
        ID3D12GraphicsCommandList *pCmdList);
    
    // Set text color (default: white)
    void setColor(const float4& color) { m_textColor = color; }

private:
    // Text vertex structure (matches shader input)
    struct TextVertex
    {
        float2 position;    // Pixel coordinates (converted to NDC in vertex shader)
        float2 texCoord;    // UV coordinates in font atlas
    };
    
    // Text parameters constant buffer structure
    struct TextParams
    {
        float4 textColor;
        float2 screenSize;
        float2 padding;
    };
    
    // Generate vertex data for all text lines
    void generateTextQuads(std::vector<TextVertex>& vertices, std::vector<uint16_t>& indices, 
                          float2 screenSize);
    
    // Update vertex buffer with new text data
    void updateVertexBuffer(const std::vector<TextVertex>& vertices, 
                           const std::vector<uint16_t>& indices, 
                           ID3D12Device* pDevice);
    
    // Create or update constant buffer
    void updateConstantBuffer(const float2& screenSize, ID3D12Device* pDevice);
    
    // Create descriptor heaps if needed
    void ensureDescriptorHeaps(ID3D12Device* pDevice);

private:
    std::shared_ptr<GPUFont> m_pFont;
    float2 m_vLeftTop;
    std::vector<std::string> m_sLines; // lines of text it's going to draw
    float4 m_textColor = float4(1.0f, 1.0f, 1.0f, 1.0f); // Default: white
    
    // Rendering resources
    Microsoft::WRL::ComPtr<ID3D12Resource> m_vertexBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_indexBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_constantBuffer;
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_descriptorHeap;  // Combined CBV/SRV heap
    
    D3D12_VERTEX_BUFFER_VIEW m_vertexBufferView = {};
    D3D12_INDEX_BUFFER_VIEW m_indexBufferView = {};
    
    uint32_t m_vertexCount = 0;
    uint32_t m_indexCount = 0;
    uint8_t* m_constantBufferData = nullptr;
};
