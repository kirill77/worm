#include "pch.h"
#include "Window.h"
#include "DirectXHelpers.h"
#include "GPUMesh.h"
#include "GPUQueue.h"

// Helper function to create a GPU buffer
Microsoft::WRL::ComPtr<ID3D12Resource> createDefaultBuffer(
    ID3D12Device* device,
    ID3D12GraphicsCommandList* cmdList,
    const void* initData,
    UINT64 byteSize,
    Microsoft::WRL::ComPtr<ID3D12Resource>& uploadBuffer)
{
    Microsoft::WRL::ComPtr<ID3D12Resource> defaultBuffer;

    // Create the actual default buffer resource
    CD3DX12_HEAP_PROPERTIES defaultHeapProperties(D3D12_HEAP_TYPE_DEFAULT);
    CD3DX12_RESOURCE_DESC bufferDesc = CD3DX12_RESOURCE_DESC::Buffer(byteSize);
    
    ThrowIfFailed(device->CreateCommittedResource(
        &defaultHeapProperties,
        D3D12_HEAP_FLAG_NONE,
        &bufferDesc,
        D3D12_RESOURCE_STATE_COMMON,
        nullptr,
        IID_PPV_ARGS(defaultBuffer.GetAddressOf())));

    // Create an upload heap to copy the data
    CD3DX12_HEAP_PROPERTIES uploadHeapProperties(D3D12_HEAP_TYPE_UPLOAD);
    
    ThrowIfFailed(device->CreateCommittedResource(
        &uploadHeapProperties,
        D3D12_HEAP_FLAG_NONE,
        &bufferDesc,
        D3D12_RESOURCE_STATE_GENERIC_READ,
        nullptr,
        IID_PPV_ARGS(uploadBuffer.GetAddressOf())));

    // Copy data to the upload heap
    D3D12_SUBRESOURCE_DATA subResourceData = {};
    subResourceData.pData = initData;
    subResourceData.RowPitch = byteSize;
    subResourceData.SlicePitch = subResourceData.RowPitch;

    // Change state and copy data
    CD3DX12_RESOURCE_BARRIER resourceBarrier = 
        CD3DX12_RESOURCE_BARRIER::Transition(
            defaultBuffer.Get(),
            D3D12_RESOURCE_STATE_COMMON,
            D3D12_RESOURCE_STATE_COPY_DEST);
    
    cmdList->ResourceBarrier(1, &resourceBarrier);
    
    // Use the non-templated version from CD3DX12.h
    UpdateSubresources(
        cmdList,
        defaultBuffer.Get(),
        uploadBuffer.Get(),
        0, 0, 1,
        &subResourceData);
    
    resourceBarrier = CD3DX12_RESOURCE_BARRIER::Transition(
        defaultBuffer.Get(),
        D3D12_RESOURCE_STATE_COPY_DEST,
        D3D12_RESOURCE_STATE_GENERIC_READ);
    
    cmdList->ResourceBarrier(1, &resourceBarrier);

    return defaultBuffer;
}

GPUMesh::GPUMesh(Microsoft::WRL::ComPtr<ID3D12Device> device)
    : m_device(device)
{
}

void GPUMesh::setGeometry(const std::vector<Vertex>& pVertices, std::vector<int3>& pTriangles, GPUQueue& gpuQueue)
{
    // Get a command list from the provided queue
    auto commandList = gpuQueue.beginRecording();
    
    // Upload buffers
    Microsoft::WRL::ComPtr<ID3D12Resource> vertexUploadBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> indexUploadBuffer;
    
    // Create vertex buffer
    const UINT vbSize = static_cast<UINT>(pVertices.size() * sizeof(Vertex));
    m_vertexBuffer = createDefaultBuffer(
        m_device.Get(),
        commandList.Get(),
        pVertices.data(),
        vbSize,
        vertexUploadBuffer);
    
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
    
    // Create index buffer
    const UINT ibSize = static_cast<UINT>(indices.size() * sizeof(UINT));
    m_indexBuffer = createDefaultBuffer(
        m_device.Get(),
        commandList.Get(),
        indices.data(),
        ibSize,
        indexUploadBuffer);
    
    // Create index buffer view
    m_indexBufferView.BufferLocation = m_indexBuffer->GetGPUVirtualAddress();
    m_indexBufferView.Format = DXGI_FORMAT_R32_UINT;
    m_indexBufferView.SizeInBytes = ibSize;
    
    // Execute command list using the provided queue
    gpuQueue.execute(commandList);
} 