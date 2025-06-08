#pragma once

#include <d3d12.h>
#include <wrl/client.h>
#include <filesystem>

// Forward declarations
namespace Microsoft { namespace WRL { template<typename> class ComPtr; } }
class GPUQueue;

struct GPUResource
{
    GPUResource(Microsoft::WRL::ComPtr<ID3D12Device> device);
    
    void loadFromFile(const std::filesystem::path& path, GPUQueue& queue);
    
    Microsoft::WRL::ComPtr<ID3D12Resource> getResource() const { return m_resource; }

private:
    Microsoft::WRL::ComPtr<ID3D12Device> m_device;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_resource;
};

