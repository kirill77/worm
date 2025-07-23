#pragma once

#include <memory>
#include <vector>
#include "geometry/vectors/vector.h"
#include "geometry/vectors/affine.h"
#include "geometry/vectors/box.h"
#include <d3d12.h>
#include <DirectXMath.h>
#include <wrl/client.h>

// Forward declarations
namespace Microsoft { namespace WRL { template<typename> class ComPtr; } }
struct ID3D12Device;
struct ID3D12Resource;

class GPUMesh
{
public:
    GPUMesh(Microsoft::WRL::ComPtr<ID3D12Device> device);
    
    struct Vertex
    {
        float3 vPos;
    };
    
    void setGeometry(const std::vector<Vertex>& pVertices, std::vector<int3>& pTriangles);
    
    D3D12_VERTEX_BUFFER_VIEW getVertexBufferView() const { return m_vertexBufferView; }
    D3D12_INDEX_BUFFER_VIEW getIndexBufferView() const { return m_indexBufferView; }
    uint32_t getIndexCount() const { return m_indexCount; }
    
    const box3& getBoundingBox() const { return m_boundingBox; }
    
    // Transform methods
    void setTransform(const affine3& transform) { m_mToParent = transform; }
    const affine3& getTransform() const { return m_mToParent; }
    DirectX::XMMATRIX getWorldMatrix() const;
    
private:
    Microsoft::WRL::ComPtr<ID3D12Device> m_device;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_vertexBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_indexBuffer;
    D3D12_VERTEX_BUFFER_VIEW m_vertexBufferView;
    D3D12_INDEX_BUFFER_VIEW m_indexBufferView;
    uint32_t m_indexCount = 0;
    box3 m_boundingBox;
    
    // Transform from mesh local space to parent space
    affine3 m_mToParent;
}; 