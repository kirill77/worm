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
    struct ButtonOrKey
    {
        uint64_t getLastChangeInputTick() const { return lastChangeInputTick; }
        
        // Access to stored Windows message data
        LPARAM getLastLParam() const { return lastLParam; }
        WPARAM getLastWParam() const { return lastWParam; }
        UINT getLastMessage() const { return lastMessage; }
        
        // Button/key state query methods
        uint16_t getRepeatCount() const;
        uint8_t getScanCode() const;
        bool isExtended() const;
        bool wasRepeated() const;
        float2 getLastClickPosition() const;

        void notifyState(UINT message, WPARAM wParam, LPARAM lParam, uint64_t inputTick);

    private:
        uint64_t lastChangeInputTick = 0; // input tick when last changed
        
        // Store full Windows message context
        LPARAM lastLParam = 0;
        WPARAM lastWParam = 0;
        UINT lastMessage = 0;
    };

    float2 getMousePosition() const;
    float getScrollWheelState() const;
    
    // Button/key state query
    bool isPressed(uint32_t keyId, bool bIgnoreRepeats = false) const;
    
    // Read-only access to ButtonOrKey struct
    const ButtonOrKey& getButtonOrKey(uint32_t buttonOrKeyId) const;
    
    // Input tick-based timing (separate from rendering frames)
    void notifyBeforeInputTick();
    uint64_t getCurrentInputTick() const;
    
    // Input update methods (previously accessed through friend relationship)
    void notifyButtonOrKeyState(UINT message, WPARAM wParam, LPARAM lParam);
    void setMousePosition(float x, float y);
    void updateScrollWheelState(float delta);

private:
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
