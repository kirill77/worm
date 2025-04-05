#pragma once

#include <memory>
#include <vector>
#include "../math/vector.h"
#include <d3d12.h>
#include <wrl/client.h>

// Forward declarations
namespace Microsoft { namespace WRL { template<typename> class ComPtr; } }
struct ID3D12Device;
struct ID3D12Resource;
class GPUQueue;

class GPUMesh
{
public:
    GPUMesh(Microsoft::WRL::ComPtr<ID3D12Device> device);
    
    struct Vertex
    {
        float3 vPos;
    };
    
    void setGeometry(const std::vector<Vertex>& pVertices, std::vector<int3>& pTriangles, GPUQueue& gpuQueue);
    
    D3D12_VERTEX_BUFFER_VIEW getVertexBufferView() const { return m_vertexBufferView; }
    D3D12_INDEX_BUFFER_VIEW getIndexBufferView() const { return m_indexBufferView; }
    uint32_t getIndexCount() const { return m_indexCount; }
    
private:
    Microsoft::WRL::ComPtr<ID3D12Device> m_device;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_vertexBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_indexBuffer;
    D3D12_VERTEX_BUFFER_VIEW m_vertexBufferView;
    D3D12_INDEX_BUFFER_VIEW m_indexBufferView;
    uint32_t m_indexCount = 0;
}; 