#pragma once

#include <memory>
#include "../math/vector.h"
#include <d3d12.h>
#include <wrl/client.h>

// Forward declarations
namespace Microsoft { namespace WRL { template<typename> class ComPtr; } }
struct ID3D12Device;
struct ID3D12CommandQueue;
struct ID3D12GraphicsCommandList;
struct ID3D12CommandAllocator;

class GPUQueue
{
public:
    GPUQueue(Microsoft::WRL::ComPtr<ID3D12Device> device);
    
    Microsoft::WRL::ComPtr<ID3D12CommandQueue> getQueue();
    Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList> beginRecording();
    bool execute(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList> pCmdList);

private:
    Microsoft::WRL::ComPtr<ID3D12Device> m_device;
    Microsoft::WRL::ComPtr<ID3D12CommandQueue> m_commandQueue;
    Microsoft::WRL::ComPtr<ID3D12CommandAllocator> m_commandAllocator;
    Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList> m_commandList;
}; 