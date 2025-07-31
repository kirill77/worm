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
    
    // Input tick-based timing (separate from rendering frames)
    void notifyBeforeInputTick();
    uint64_t getCurrentInputTick() const;
    uint64_t getLastChangeInputTick(uint32_t buttonOrKeyId) const;
    uint64_t getInputTicksSinceLastChange(uint32_t buttonOrKeyId) const;
    bool wasChangedInInputTick(uint32_t buttonOrKeyId, uint64_t inputTick) const;
    bool wasChangedInLastNInputTicks(uint32_t buttonOrKeyId, uint64_t tickCount) const;
    
    // Enhanced query methods for rich input information
    uint16_t getKeyRepeatCount(uint32_t keyId) const;
    uint8_t getKeyScanCode(uint32_t keyId) const;
    bool isExtendedKey(uint32_t keyId) const;
    bool wasKeyRepeated(uint32_t keyId) const;
    float2 getLastClickPosition(uint32_t mouseButton) const;
    
    // Input update methods (previously accessed through friend relationship)
    void notifyButtonOrKeyState(UINT message, WPARAM wParam, LPARAM lParam);
    void setMousePosition(float x, float y);
    void updateScrollWheelState(float delta);

private:
    struct ButtonOrKey
    {
        void notifyPressed(uint64_t inputTick);
        void notifyReleased(uint64_t inputTick);
        void notifyState(UINT message, WPARAM wParam, LPARAM lParam, uint64_t inputTick);
        uint32_t getPressCount() const;
        uint32_t getReleaseCount() const;
        uint64_t getLastChangeInputTick() const { return lastChangeInputTick; }
        
        // Access to stored Windows message data
        WPARAM getLastWParam() const { return lastWParam; }
        LPARAM getLastLParam() const { return lastLParam; }
        UINT getLastMessage() const { return lastMessage; }
        
    private:
        uint32_t pressCount = 0; // how many times was the control pressed
        uint32_t releaseCount = 0; // how many times was the control released
        uint64_t lastChangeInputTick = 0; // input tick when last changed
        
        // Store full Windows message context
        WPARAM lastWParam = 0;
        LPARAM lastLParam = 0;
        UINT lastMessage = 0;
    };
    std::unordered_map<uint32_t, ButtonOrKey> m_buttonsAndKeys;
    float2 m_mousePosition = float2(0.0f, 0.0f);
    float m_scrollWheelState = 0.0f;
    uint64_t m_currentInputTick = 0;
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
