#pragma once

#include <memory>
#include <d3d12.h>
#include <dxgi1_6.h>
#include <wrl/client.h>
#include "GPUQueue.h"

struct SwapChain : public std::enable_shared_from_this<SwapChain>
{
    SwapChain(ID3D12Device *pDevice, HWND hWnd);

    IDXGISwapChain4* getSwapChain();
    ID3D12CommandQueue* getCommandQueue();
    GPUQueue* getGPUQueue() { return m_pGPUQueue.get(); }
    void notifyWindowResized();

    // access to current back-buffer resource and handles
    ID3D12Resource* getBBColor();
    D3D12_CPU_DESCRIPTOR_HANDLE getBBColorCPUHandle();
    D3D12_GPU_DESCRIPTOR_HANDLE getBBColorGPUHandle();
    ID3D12Resource* getBBDepth();
    D3D12_CPU_DESCRIPTOR_HANDLE getBBDepthCPUHandle();
    D3D12_GPU_DESCRIPTOR_HANDLE getBBDepthGPUHandle();

private:
    void createBackBufferResources();
    void createDepthBuffer();
    void releaseBackBufferResources();

    Microsoft::WRL::ComPtr<ID3D12Device> m_pDevice;
    std::shared_ptr<GPUQueue> m_pGPUQueue;
    Microsoft::WRL::ComPtr<IDXGISwapChain4> m_pSwapChain;
    HWND m_hWnd;
    
    // Descriptor heaps
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_pRTVHeap;
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_pDSVHeap;
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_pSRVHeap;
    
    // Back buffer resources
    static const UINT m_backBufferCount = 2;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_pBackBuffers[m_backBufferCount];
    Microsoft::WRL::ComPtr<ID3D12Resource> m_pDepthBuffer;
    
    UINT m_rtvDescriptorSize;
    UINT m_dsvDescriptorSize;
    UINT m_srvDescriptorSize;
};

