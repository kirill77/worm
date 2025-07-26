#pragma once

#include <memory>
#include <vector>
#include <DirectXMath.h>
#include <d3d12.h>
#include <wrl/client.h>
#include "GPUCamera.h"
#include "GPUFont.h"
#include "IObjectVis.h"
#include "geometry/vectors/box.h"

// Forward declarations
namespace Microsoft { namespace WRL { template<typename> class ComPtr; } }
struct ID3D12RootSignature;
struct ID3D12PipelineState;
struct ID3D12Resource;
struct ID3D12DescriptorHeap;
struct Window;
struct GPUStats;

struct GPUWorld
{
public:
    GPUWorld(std::shared_ptr<Window> pWindow, GPUQueue* pGpuQueue);
    ~GPUWorld();

    // Object management
    void addObject(std::weak_ptr<IObjectVis> pObject);
    
    // Camera management
    std::shared_ptr<GPUCamera> getCamera();
    void setCamera(std::shared_ptr<GPUCamera> camera);
    
    // Font management
    std::shared_ptr<GPUFont> getFont();
    
    // Root signature access (shared between mesh and text rendering)
    ID3D12RootSignature* getSharedRootSignature() const { return m_pRootSignature.Get(); }
    
    // Rendering - returns bounding box of all visualized objects
    box3 render(SwapChain* pSwapChain, ID3D12GraphicsCommandList* pCmdList);

private:
    void initializeRenderResources();
    
    struct TransformBuffer {
        DirectX::XMMATRIX View;
        DirectX::XMMATRIX Projection;
    };
    
    std::shared_ptr<Window> m_pWindow;
    // when weak_ptr becomes null - it's automatically removed from the list
    std::vector<std::weak_ptr<IObjectVis>> m_pObjects;
    std::shared_ptr<GPUCamera> m_pCamera;
    std::shared_ptr<GPUFont> m_pFont;
    
    // DirectX rendering resources
    Microsoft::WRL::ComPtr<ID3D12RootSignature> m_pRootSignature;
    Microsoft::WRL::ComPtr<ID3D12PipelineState> m_pPipelineState;

    Microsoft::WRL::ComPtr<ID3D12Resource> m_pTransformRes;
    TransformBuffer m_transformBufferMatrix;
    UINT8* m_pTransformData = nullptr;

    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_pCBVHeap;
}; 