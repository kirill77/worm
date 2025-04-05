#pragma once

#include <memory>
#include <vector>
#include <DirectXMath.h>
#include <d3d12.h>
#include <wrl/client.h>
#include "GPUMesh.h"
#include "GPUCamera.h"

// Forward declarations
namespace Microsoft { namespace WRL { template<typename> class ComPtr; } }
struct ID3D12RootSignature;
struct ID3D12PipelineState;
struct ID3D12Resource;
struct ID3D12DescriptorHeap;
class Window;
class GPUStats;

class GPUWorld
{
public:
    GPUWorld(std::shared_ptr<Window> pWindow);
    ~GPUWorld();

    // Mesh management
    std::shared_ptr<GPUMesh> createMesh();
    void addMesh(std::shared_ptr<GPUMesh> mesh);
    void removeMesh(std::shared_ptr<GPUMesh> mesh);
    
    // Camera management
    std::shared_ptr<GPUCamera> getCamera();
    void setCamera(std::shared_ptr<GPUCamera> camera);
    
    // Rendering
    void drawMeshesIntoWindow(GPUStats *pStats = nullptr);

private:
    void initializeRenderResources();
    
    struct TransformBuffer {
        DirectX::XMMATRIX World;
        DirectX::XMMATRIX View;
        DirectX::XMMATRIX Projection;
    };
    
    std::shared_ptr<Window> m_pWindow;
    std::vector<std::shared_ptr<GPUMesh>> m_pMeshes;
    std::shared_ptr<GPUCamera> m_pCamera;
    
    // DirectX rendering resources
    Microsoft::WRL::ComPtr<ID3D12RootSignature> m_rootSignature;
    Microsoft::WRL::ComPtr<ID3D12PipelineState> m_pipelineState;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_transformBufferResource;
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_cbvHeap;
    TransformBuffer m_transformBufferMatrix;
    UINT8* m_transformBufferData = nullptr;
}; 