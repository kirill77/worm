#include "pch.h"
#include "Window.h"
#include "DirectXHelpers.h"
#include <stdexcept>

#pragma comment(lib, "d3d12.lib")
#pragma comment(lib, "dxgi.lib")
#pragma comment(lib, "d3dcompiler.lib")

// Constants for window creation
constexpr int DEFAULT_WIDTH = 1280;
constexpr int DEFAULT_HEIGHT = 720;

// Global Window class pointer for the Win32 window procedure
Window* g_pWindow = nullptr;

// ButtonOrKey implementation
void UIState::ButtonOrKey::notifyPressed() {
    pressCount++;
    lastChangeTS = std::time(nullptr);
}

void UIState::ButtonOrKey::notifyReleased() {
    releaseCount++;
    lastChangeTS = std::time(nullptr);
}

void UIState::ButtonOrKey::notifyState(UINT message, WPARAM wParam, LPARAM lParam) {
    // Update counters based on message type
    if (message == WM_KEYDOWN || message == WM_LBUTTONDOWN || message == WM_RBUTTONDOWN || message == WM_MBUTTONDOWN) {
        pressCount++;
    } else if (message == WM_KEYUP || message == WM_LBUTTONUP || message == WM_RBUTTONUP || message == WM_MBUTTONUP) {
        releaseCount++;
    }
    
    // Store full Windows message context
    lastWParam = wParam;
    lastLParam = lParam;
    lastMessage = message;
    lastChangeTS = std::time(nullptr);
}

uint32_t UIState::ButtonOrKey::getPressCount() const {
    return pressCount;
}

uint32_t UIState::ButtonOrKey::getReleaseCount() const {
    return releaseCount;
}

// UIState implementation
bool UIState::isButtonOrKeyPressed(uint32_t buttonOrKeyId) const {
    auto it = m_buttonsAndKeys.find(buttonOrKeyId);
    if (it == m_buttonsAndKeys.end()) return false;
    return (GetAsyncKeyState(buttonOrKeyId) & 0x8000) != 0;
}

uint32_t UIState::getButtonOrKeyPressCount(uint32_t buttonOrKeyId) const {
    auto it = m_buttonsAndKeys.find(buttonOrKeyId);
    if (it != m_buttonsAndKeys.end()) {
        return it->second.getPressCount();
    }
    return 0;
}

float2 UIState::getMousePosition() const {
    return m_mousePosition;
}

float UIState::getScrollWheelState() const {
    return m_scrollWheelState;
}

void UIState::notifyButtonOrKeyState(UINT message, WPARAM wParam, LPARAM lParam) {
    uint32_t keyId;
    
    // Determine key/button ID from message type and wParam
    switch (message) {
        case WM_KEYDOWN:
        case WM_KEYUP:
            keyId = static_cast<uint32_t>(wParam);
            break;
        case WM_LBUTTONDOWN:
        case WM_LBUTTONUP:
            keyId = VK_LBUTTON;
            break;
        case WM_RBUTTONDOWN:
        case WM_RBUTTONUP:
            keyId = VK_RBUTTON;
            break;
        case WM_MBUTTONDOWN:
        case WM_MBUTTONUP:
            keyId = VK_MBUTTON;
            break;
        default:
            return; // Unknown message type
    }
    
    m_buttonsAndKeys[keyId].notifyState(message, wParam, lParam);
}

void UIState::setMousePosition(float x, float y) {
    m_mousePosition.x = x;
    m_mousePosition.y = y;
}

void UIState::updateScrollWheelState(float delta) {
    m_scrollWheelState += delta;
}

// Enhanced query methods for rich input information
uint16_t UIState::getKeyRepeatCount(uint32_t keyId) const {
    auto it = m_buttonsAndKeys.find(keyId);
    if (it != m_buttonsAndKeys.end() && it->second.getLastMessage() == WM_KEYDOWN) {
        return static_cast<uint16_t>(it->second.getLastLParam() & 0xFFFF);
    }
    return 0;
}

uint8_t UIState::getKeyScanCode(uint32_t keyId) const {
    auto it = m_buttonsAndKeys.find(keyId);
    if (it != m_buttonsAndKeys.end()) {
        return static_cast<uint8_t>((it->second.getLastLParam() >> 16) & 0xFF);
    }
    return 0;
}

bool UIState::isExtendedKey(uint32_t keyId) const {
    auto it = m_buttonsAndKeys.find(keyId);
    if (it != m_buttonsAndKeys.end()) {
        return (it->second.getLastLParam() & (1 << 24)) != 0;
    }
    return false;
}

bool UIState::wasKeyRepeated(uint32_t keyId) const {
    auto it = m_buttonsAndKeys.find(keyId);
    if (it != m_buttonsAndKeys.end()) {
        return (it->second.getLastLParam() & (1 << 30)) != 0; // Previous key state
    }
    return false;
}

float2 UIState::getLastClickPosition(uint32_t mouseButton) const {
    auto it = m_buttonsAndKeys.find(mouseButton);
    if (it != m_buttonsAndKeys.end()) {
        UINT msg = it->second.getLastMessage();
        if (msg == WM_LBUTTONDOWN || msg == WM_RBUTTONDOWN || msg == WM_MBUTTONDOWN) {
            LPARAM lParam = it->second.getLastLParam();
            return float2(static_cast<float>(GET_X_LPARAM(lParam)),
                         static_cast<float>(GET_Y_LPARAM(lParam)));
        }
    }
    return float2(0.0f, 0.0f);
}

// Window class implementation
Window::Window() 
    : m_hwnd(nullptr)
    , m_width(DEFAULT_WIDTH)
    , m_height(DEFAULT_HEIGHT)
    , m_uiState(std::make_unique<UIState>())
{
    g_pWindow = this;
}

Window::~Window() {
    // Release DirectX resources
    m_pSwapChain.reset();
    m_device.Reset();
    
    // Unregister window class
    if (m_hwnd) {
        DestroyWindow(m_hwnd);
        m_hwnd = nullptr;
    }
    
    g_pWindow = nullptr;
}

bool Window::createWindowDevicAndSwapChain(const std::string& sName) {
    // Disable Windows DPI scaling - we'll handle scaling ourselves
    // Try the newer API first (Windows 10 version 1703+)
    typedef BOOL(WINAPI* SetProcessDpiAwarenessContextProc)(DPI_AWARENESS_CONTEXT);
    HMODULE user32 = GetModuleHandle(L"user32.dll");
    if (user32) {
        SetProcessDpiAwarenessContextProc setProcessDpiAwarenessContext = 
            (SetProcessDpiAwarenessContextProc)GetProcAddress(user32, "SetProcessDpiAwarenessContext");
        if (setProcessDpiAwarenessContext) {
            // Tell Windows we are per-monitor DPI aware but will handle scaling ourselves
            setProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);
        } else {
            // Fallback for older Windows versions
            SetProcessDPIAware();
        }
    }

    // Get desktop resolution and set window size to match
    m_width = GetSystemMetrics(SM_CXSCREEN);
    m_height = GetSystemMetrics(SM_CYSCREEN);

    // Register window class
    WNDCLASSEX windowClass = {};
    windowClass.cbSize = sizeof(WNDCLASSEX);
    windowClass.style = CS_HREDRAW | CS_VREDRAW;
    windowClass.lpfnWndProc = Window::WindowProc;
    windowClass.hInstance = GetModuleHandle(NULL);
    windowClass.hCursor = LoadCursor(NULL, IDC_ARROW);
    windowClass.lpszClassName = L"VisualizationWindowClass";
    RegisterClassEx(&windowClass);

    // Create window
    RECT windowRect = { 0, 0, static_cast<LONG>(m_width), static_cast<LONG>(m_height) };
    AdjustWindowRect(&windowRect, WS_POPUP, FALSE);

    // Convert string name to wstring
    std::wstring wideName(sName.begin(), sName.end());

    m_hwnd = CreateWindow(
        windowClass.lpszClassName,
        wideName.c_str(),
        WS_POPUP,
        CW_USEDEFAULT,
        CW_USEDEFAULT,
        windowRect.right - windowRect.left,
        windowRect.bottom - windowRect.top,
        nullptr, // No parent window
        nullptr, // No menus
        windowClass.hInstance,
        this);

    if (!m_hwnd) {
        return false;
    }

    // Show the window
    ShowWindow(m_hwnd, SW_SHOW);

    // Initialize DirectX
    if (!initDirectX()) {
        return false;
    }

    return true;
}

const UIState& Window::getCurrentUIState() {
    return *m_uiState;
}


void Window::processMessages() {
    MSG msg = {};
    
    // Process all pending Windows messages
    while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
}

HWND Window::getWindowHandle() const {
    return m_hwnd;
}

bool Window::initDirectX() {
    UINT dxgiFactoryFlags = 0;

#ifdef _DEBUG
    // Enable debug layer
    Microsoft::WRL::ComPtr<ID3D12Debug> debugController;
    if (SUCCEEDED(D3D12GetDebugInterface(IID_PPV_ARGS(&debugController)))) {
        debugController->EnableDebugLayer();
        dxgiFactoryFlags |= DXGI_CREATE_FACTORY_DEBUG;
    }
#endif

    // Create DXGI factory
    Microsoft::WRL::ComPtr<IDXGIFactory6> dxgiFactory;
    ThrowIfFailed(CreateDXGIFactory2(dxgiFactoryFlags, IID_PPV_ARGS(&dxgiFactory)));

    // Create D3D12 device
    Microsoft::WRL::ComPtr<IDXGIAdapter1> hardwareAdapter;
    for (UINT adapterIndex = 0; DXGI_ERROR_NOT_FOUND != dxgiFactory->EnumAdapters1(adapterIndex, &hardwareAdapter); ++adapterIndex) {
        DXGI_ADAPTER_DESC1 adapterDesc;
        hardwareAdapter->GetDesc1(&adapterDesc);

        // Skip software adapters
        if (adapterDesc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE) {
            continue;
        }

        // Try to create the device
        if (SUCCEEDED(D3D12CreateDevice(hardwareAdapter.Get(), D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&m_device)))) {
            break;
        }
    }

    if (!m_device) {
        // If no hardware adapter found, use WARP adapter
        ThrowIfFailed(dxgiFactory->EnumWarpAdapter(IID_PPV_ARGS(&hardwareAdapter)));
        ThrowIfFailed(D3D12CreateDevice(hardwareAdapter.Get(), D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&m_device)));
    }

    // Create SwapChain object (which creates and manages the command queue)
    m_pSwapChain = std::make_shared<SwapChain>(m_device.Get(), m_hwnd);

    // Disable Alt+Enter fullscreen toggle
    ThrowIfFailed(dxgiFactory->MakeWindowAssociation(m_hwnd, DXGI_MWA_NO_ALT_ENTER));

    return true;
}

// Helper function to handle mouse/keyboard input
void Window::handleInput(UINT message, WPARAM wParam, LPARAM lParam) {
    switch (message) {
    case WM_KEYDOWN:
    case WM_KEYUP:
    case WM_LBUTTONDOWN:
    case WM_LBUTTONUP:
    case WM_RBUTTONDOWN:
    case WM_RBUTTONUP:
    case WM_MBUTTONDOWN:
    case WM_MBUTTONUP:
        // Use the new method that captures full Windows message information
        m_uiState->notifyButtonOrKeyState(message, wParam, lParam);
        break;
    case WM_MOUSEMOVE:
        m_uiState->setMousePosition(static_cast<float>(GET_X_LPARAM(lParam)), 
                                     static_cast<float>(GET_Y_LPARAM(lParam)));
        break;
    case WM_MOUSEWHEEL:
        m_uiState->updateScrollWheelState(static_cast<float>(GET_WHEEL_DELTA_WPARAM(wParam)) / WHEEL_DELTA);
        break;
    }
}

// Window procedure
LRESULT CALLBACK Window::WindowProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    Window* window = nullptr;
    if (message == WM_NCCREATE)
    {
        CREATESTRUCT* createStruct = reinterpret_cast<CREATESTRUCT*>(lParam);
        window = reinterpret_cast<Window*>(createStruct->lpCreateParams);
        SetWindowLongPtr(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(window));
    }
    else
    {
        window = reinterpret_cast<Window*>(GetWindowLongPtr(hwnd, GWLP_USERDATA));
    }

    if (window)
    {
        switch (message)
        {
        case WM_DESTROY:
            window->m_shouldExit = true;
            return 0;
        case WM_SIZE:
        {
            // Get the new window dimensions
            window->m_width = LOWORD(lParam);
            window->m_height = HIWORD(lParam);

            // Notify GPUWorld about the resize
            window->onWindowResize(LOWORD(lParam), HIWORD(lParam));
            return 0;
        }
        case WM_KEYDOWN:
        case WM_KEYUP:
        case WM_LBUTTONDOWN:
        case WM_LBUTTONUP:
        case WM_RBUTTONDOWN:
        case WM_RBUTTONUP:
        case WM_MOUSEMOVE:
            window->handleInput(message, wParam, lParam);
            return 0;
        }
    }

    return DefWindowProc(hwnd, message, wParam, lParam);
}

void Window::onWindowResize(UINT width, UINT height)
{
    if (width == 0 || height == 0)
        return;

    if (m_pSwapChain)
    {
        // Wait for GPU to complete all operations
        GPUQueue* gpuQueue = m_pSwapChain->getGPUQueue();
        if (gpuQueue)
        {
            gpuQueue->flush();
        }

        // Use SwapChain's resize functionality
        m_pSwapChain->notifyWindowResized();
    }

    // Update window dimensions
    m_width = width;
    m_height = height;
}
