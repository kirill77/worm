#pragma once

#include <memory>
#include <cstdarg>
#include <vector>
#include <ctime>
#include "GPUFont.h"
#include "geometry/vectors/vector.h"
#include <d3d12.h>
#include <wrl/client.h>

// Forward declarations
struct GPUQueue;
struct SwapChain;

// Line data structure for text lines
struct Line
{
public:
    Line() : m_createTS(std::time(nullptr)), m_color(1.0f, 1.0f, 1.0f, 1.0f) {}
    
    // Set the text content using printf-style formatting
    int printf(const char* format, ...);
    
    // Getters
    const std::string& getText() const { return m_string; }
    bool isEmpty() const { return m_string.empty(); }
    const float4& getColor() const { return m_color; }
    uint32_t getLifeTimeSec() const { return m_lifeTimeSec; }
    std::time_t getCreateTime() const { return m_createTS; }
    
    // Setters with validation
    void setColor(const float4& color);
    void setLifeTime(uint32_t lifeTimeSec) { m_lifeTimeSec = lifeTimeSec; }

private:
    std::string m_string;
    std::time_t m_createTS; // when the line first appeared
    uint32_t m_lifeTimeSec = 0; // how many seconds to keep it on the screen, 0 if forever
    float4 m_color; // default - white
};

struct GPUText
{
public:
    GPUText(std::shared_ptr<GPUFont> pFont);

    void setLeftTop(const float2& vLeftTop);
    
    // Create a new line and return a shared pointer to it
    std::shared_ptr<Line> createLine();
    
    // Check if a line has expired based on its lifetime
    bool isExpired(const std::shared_ptr<Line>& line) const;
    
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
    std::vector<std::shared_ptr<Line>> m_pLines; // lines of text it's going to draw
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
