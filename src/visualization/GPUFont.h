#pragma once

#include <memory>
#include <unordered_map>
#include "GPUResource.h"
#include "math/vector.h"
#include <d3d12.h>
#include <wrl/client.h>

// Forward declarations
struct GPUQueue;
struct ID3D12RootSignature;
struct ID3D12PipelineState;

struct GlyphInfo {
    float2 texCoords[2];    // UV coordinates in atlas (min, max)
    float2 size;            // Glyph dimensions in pixels
    float2 bearing;         // Offset from baseline
    float advance;          // Horizontal advance to next character
};

struct GPUFont : public std::enable_shared_from_this<GPUFont>
{
    GPUFont(uint32_t fontSize, GPUQueue* pQueue);

    GPUResource* getResource() const { return m_pFont.get(); }
    const GlyphInfo* getGlyphInfo(char character) const;
    float getLineHeight() const { return m_lineHeight; }
    
    // Get the text rendering PSO (creates it if needed)
    ID3D12PipelineState* getTextPSO(ID3D12RootSignature* pRootSignature);

protected:
    // Create the pipeline state object for text rendering
    void createPSO(ID3D12RootSignature* pRootSignature);

private:
    std::shared_ptr<GPUResource> m_pFont;
    std::unordered_map<char, GlyphInfo> m_glyphMap;
    float m_fontSize;
    float m_lineHeight;
    
    // Text rendering PSO
    Microsoft::WRL::ComPtr<ID3D12PipelineState> m_textPSO;
    Microsoft::WRL::ComPtr<ID3D12Device> m_device;
};

