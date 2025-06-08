#include "pch.h"
#include "Window.h"
#include "DirectXHelpers.h"
#include "GPUMesh.h"
#include <stdexcept>

// Helper function to create or update a GPU buffer in CPU-visible memory
static Microsoft::WRL::ComPtr<ID3D12Resource> createOrUpdateUploadBuffer(
    ID3D12Device* device,
    const void* initData,
    UINT64 byteSize,
    Microsoft::WRL::ComPtr<ID3D12Resource>& existingBuffer)
{
    Microsoft::WRL::ComPtr<ID3D12Resource> uploadBuffer;

    // If buffer already exists, validate size and update it
    if (existingBuffer)
    {
        D3D12_RESOURCE_DESC existingDesc = existingBuffer->GetDesc();
        if (existingDesc.Width != byteSize)
        {
            throw std::runtime_error("Buffer size changed in createOrUpdateUploadBuffer. Expected: " 
                + std::to_string(existingDesc.Width) + ", Got: " + std::to_string(byteSize));
        }
        
        // Use existing buffer
        uploadBuffer = existingBuffer;
    }
    else
    {
        // Create buffer in upload heap (CPU-visible)
        CD3DX12_HEAP_PROPERTIES uploadHeapProperties(D3D12_HEAP_TYPE_UPLOAD);
        CD3DX12_RESOURCE_DESC bufferDesc = CD3DX12_RESOURCE_DESC::Buffer(byteSize);
        
        ThrowIfFailed(device->CreateCommittedResource(
            &uploadHeapProperties,
            D3D12_HEAP_FLAG_NONE,
            &bufferDesc,
            D3D12_RESOURCE_STATE_GENERIC_READ,
            nullptr,
            IID_PPV_ARGS(uploadBuffer.GetAddressOf())));
    }

    // Map and copy data directly
    void* mappedData = nullptr;
    ThrowIfFailed(uploadBuffer->Map(0, nullptr, &mappedData));
    memcpy(mappedData, initData, static_cast<size_t>(byteSize));
    uploadBuffer->Unmap(0, nullptr);

    return uploadBuffer;
}



GPUMesh::GPUMesh(Microsoft::WRL::ComPtr<ID3D12Device> device)
    : m_device(device)
{
}

void GPUMesh::setGeometry(const std::vector<Vertex>& pVertices, std::vector<int3>& pTriangles)
{
    // Create vertex buffer using CPU-visible upload heap
    const UINT vbSize = static_cast<UINT>(pVertices.size() * sizeof(Vertex));
    m_vertexBuffer = createOrUpdateUploadBuffer(
        m_device.Get(),
        pVertices.data(),
        vbSize,
        m_vertexBuffer);
    
    // Create vertex buffer view
    m_vertexBufferView.BufferLocation = m_vertexBuffer->GetGPUVirtualAddress();
    m_vertexBufferView.StrideInBytes = sizeof(Vertex);
    m_vertexBufferView.SizeInBytes = vbSize;
    
    // Convert int3 indices to UINT array
    std::vector<UINT> indices;
    for (const auto& tri : pTriangles)
    {
        indices.push_back(static_cast<UINT>(tri.x));
        indices.push_back(static_cast<UINT>(tri.y));
        indices.push_back(static_cast<UINT>(tri.z));
    }
    
    m_indexCount = static_cast<UINT>(indices.size());
    
    // Create index buffer using CPU-visible upload heap
    const UINT ibSize = static_cast<UINT>(indices.size() * sizeof(UINT));
    m_indexBuffer = createOrUpdateUploadBuffer(
        m_device.Get(),
        indices.data(),
        ibSize,
        m_indexBuffer);
    
    // Create index buffer view
    m_indexBufferView.BufferLocation = m_indexBuffer->GetGPUVirtualAddress();
    m_indexBufferView.Format = DXGI_FORMAT_R32_UINT;
    m_indexBufferView.SizeInBytes = ibSize;
} 