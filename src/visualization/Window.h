#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include "../math/vector.h"
#include "GPUQueue.h"
#include <d3d12.h>
#include <dxgi1_6.h>
#include <wrl/client.h>

// Forward declarations
namespace Microsoft { namespace WRL { template<typename> class ComPtr; } }
struct ID3D12Device;
struct IDXGISwapChain4;
typedef struct HWND__* HWND;

class UIState
{
public:
    friend class Window;
    
    uint32_t getButtonOrKeyPressCount(uint32_t buttonOrKeyId) const;
    float2 getMousePosition();
    float getScrollWheelState();

private:
    std::unordered_map<uint32_t, uint32_t> m_buttonKeyPressCount;
    float2 m_mousePosition = float2(0.0f, 0.0f);
    float m_scrollWheelState = 0.0f;
};

class Window
{
public:
    Window();
    ~Window();

    bool createWindowDevicAndSwapChain(const std::string &sName);
    const UIState &getCurrentUIState();

    Microsoft::WRL::ComPtr<ID3D12Device> getDevice();
    Microsoft::WRL::ComPtr<IDXGISwapChain4> getSwapChain();
    std::shared_ptr<GPUQueue> createOrGetGPUQueue();

    void processMessages();
    HWND getWindowHandle() const;
    void handleInput(UINT message, WPARAM wParam, LPARAM lParam);
    bool shouldExit() const { return m_shouldExit; }

private:
    static LRESULT CALLBACK WindowProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam);

    bool initDirectX();
    
    HWND m_hwnd;
    uint32_t m_width;
    uint32_t m_height;
    std::unique_ptr<UIState> m_uiState;
    
    Microsoft::WRL::ComPtr<ID3D12Device> m_device;
    Microsoft::WRL::ComPtr<IDXGISwapChain4> m_swapChain;
    std::shared_ptr<GPUQueue> m_gpuQueue;
    bool m_shouldExit = false;
};
