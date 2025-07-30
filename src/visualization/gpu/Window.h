#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include "geometry/vectors/vector.h"
#include "GPUQueue.h"
#include "SwapChain.h"
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
    bool isButtonOrKeyPressed(uint32_t buttonOrKeyId) const;
    uint32_t getButtonOrKeyPressCount(uint32_t buttonOrKeyId) const;
    float2 getMousePosition() const;
    float getScrollWheelState() const;
    
    // Input update methods (previously accessed through friend relationship)
    void notifyButtonOrKeyPressed(uint32_t buttonOrKeyId);
    void notifyButtonOrKeyReleased(uint32_t buttonOrKeyId);
    void setMousePosition(float x, float y);
    void updateScrollWheelState(float delta);

private:
    struct ButtonOrKey
    {
        void notifyPressed();
        void notifyReleased();
        uint32_t getPressCount() const;
        uint32_t getReleaseCount() const;
    private:
        uint32_t pressCount = 0; // how many times was the control pressed
        uint32_t releaseCount = 0; // how many times was the control released
        std::time_t lastChangeTS = 0; // used to remove old entries
    };
    std::unordered_map<uint32_t, ButtonOrKey> m_buttonsAndKeys;
    float2 m_mousePosition = float2(0.0f, 0.0f);
    float m_scrollWheelState = 0.0f;
};

struct Window
{
public:
    Window();
    ~Window();

    bool createWindowDevicAndSwapChain(const std::string &sName);
    const UIState &getCurrentUIState();

    ID3D12Device* getDevice() { return m_device.Get(); }
    SwapChain* getSwapChain() { return m_pSwapChain.get(); }

    void processMessages();
    HWND getWindowHandle() const;
    void handleInput(UINT message, WPARAM wParam, LPARAM lParam);
    bool shouldExit() const { return m_shouldExit; }
    void onWindowResize(UINT width, UINT height);

private:
    static LRESULT CALLBACK WindowProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam);

    bool initDirectX();
    
    HWND m_hwnd;
    uint32_t m_width;
    uint32_t m_height;
    std::unique_ptr<UIState> m_uiState;
    
    Microsoft::WRL::ComPtr<ID3D12Device> m_device;

    std::shared_ptr<SwapChain> m_pSwapChain;

    bool m_shouldExit = false;
};
