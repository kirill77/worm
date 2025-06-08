#include "pch.h"
#include "GPUResource.h"
#include "GPUQueue.h"
#include "DirectXHelpers.h"
#define STB_IMAGE_IMPLEMENTATION
#include "external/stb/stb_image.h"

GPUResource::GPUResource(Microsoft::WRL::ComPtr<ID3D12Device> device)
    : m_device(device)
{
}

void GPUResource::loadFromFile(const std::filesystem::path& path, GPUQueue& queue)
{
    // Load image using STB Image
    int width, height, channels;
    unsigned char* imageData = stbi_load(path.string().c_str(), &width, &height, &channels, STBI_rgb_alpha);
    if (!imageData)
    {
        throw std::runtime_error("Failed to load image: " + path.string());
    }

    // Calculate image size in bytes (RGBA format)
    const UINT imageSize = width * height * 4;

    // Create the texture resource
    D3D12_RESOURCE_DESC textureDesc = {};
    textureDesc.MipLevels = 1;
    textureDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    textureDesc.Width = width;
    textureDesc.Height = height;
    textureDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
    textureDesc.DepthOrArraySize = 1;
    textureDesc.SampleDesc.Count = 1;
    textureDesc.SampleDesc.Quality = 0;
    textureDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;

    CD3DX12_HEAP_PROPERTIES defaultHeapProperties(D3D12_HEAP_TYPE_DEFAULT);
    ThrowIfFailed(m_device->CreateCommittedResource(
        &defaultHeapProperties,
        D3D12_HEAP_FLAG_NONE,
        &textureDesc,
        D3D12_RESOURCE_STATE_COPY_DEST,
        nullptr,
        IID_PPV_ARGS(&m_resource)));

    // Create upload buffer - calculate size manually
    const UINT64 uploadBufferSize = imageSize;
    Microsoft::WRL::ComPtr<ID3D12Resource> uploadBuffer;

    CD3DX12_HEAP_PROPERTIES uploadHeapProperties(D3D12_HEAP_TYPE_UPLOAD);
    CD3DX12_RESOURCE_DESC uploadBufferDesc = CD3DX12_RESOURCE_DESC::Buffer(uploadBufferSize);
    
    ThrowIfFailed(m_device->CreateCommittedResource(
        &uploadHeapProperties,
        D3D12_HEAP_FLAG_NONE,
        &uploadBufferDesc,
        D3D12_RESOURCE_STATE_GENERIC_READ,
        nullptr,
        IID_PPV_ARGS(&uploadBuffer)));

    // Get command list from queue
    auto commandList = queue.beginRecording();

    // Map the upload buffer and copy image data
    UINT8* pData;
    CD3DX12_RANGE readRange(0, 0); // We do not intend to read from this resource on the CPU
    ThrowIfFailed(uploadBuffer->Map(0, &readRange, reinterpret_cast<void**>(&pData)));
    memcpy(pData, imageData, imageSize);
    uploadBuffer->Unmap(0, nullptr);

    // Copy from upload buffer to texture
    D3D12_TEXTURE_COPY_LOCATION srcLocation = {};
    srcLocation.pResource = uploadBuffer.Get();
    srcLocation.Type = D3D12_TEXTURE_COPY_TYPE_PLACED_FOOTPRINT;
    srcLocation.PlacedFootprint.Offset = 0;
    srcLocation.PlacedFootprint.Footprint.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    srcLocation.PlacedFootprint.Footprint.Width = width;
    srcLocation.PlacedFootprint.Footprint.Height = height;
    srcLocation.PlacedFootprint.Footprint.Depth = 1;
    srcLocation.PlacedFootprint.Footprint.RowPitch = width * 4;

    D3D12_TEXTURE_COPY_LOCATION dstLocation = {};
    dstLocation.pResource = m_resource.Get();
    dstLocation.Type = D3D12_TEXTURE_COPY_TYPE_SUBRESOURCE_INDEX;
    dstLocation.SubresourceIndex = 0;

    commandList->CopyTextureRegion(&dstLocation, 0, 0, 0, &srcLocation, nullptr);

    // Transition resource to shader resource state
    CD3DX12_RESOURCE_BARRIER barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        m_resource.Get(),
        D3D12_RESOURCE_STATE_COPY_DEST,
        D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE);
    
    commandList->ResourceBarrier(1, &barrier);

    // Execute command list
    queue.execute(commandList);

    // Free the loaded image data
    stbi_image_free(imageData);
}
