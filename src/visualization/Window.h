#pragma once

#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include "../math/vector.h"

// Forward declarations
namespace Microsoft { namespace WRL { template<typename> class ComPtr; } }
struct ID3D12Device;
struct ID3D12CommandQueue;
struct ID3D12GraphicsCommandList;
struct IDXGISwapChain4;
struct ID3D12RootSignature;
struct ID3D12PipelineState;
struct ID3D12Resource;
struct ID3D12DescriptorHeap;
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

class GPUMesh
{
public:
    GPUMesh(Microsoft::WRL::ComPtr<ID3D12Device> device);
    
    struct Vertex
    {
        float3 vPos;
    };
    
    void setGeometry(const std::vector<Vertex>& pVertices, std::vector<int3>& pTriangles);
    
    D3D12_VERTEX_BUFFER_VIEW getVertexBufferView() const { return m_vertexBufferView; }
    D3D12_INDEX_BUFFER_VIEW getIndexBufferView() const { return m_indexBufferView; }
    uint32_t getIndexCount() const { return m_indexCount; }
    
private:
    Microsoft::WRL::ComPtr<ID3D12Device> m_device;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_vertexBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_indexBuffer;
    D3D12_VERTEX_BUFFER_VIEW m_vertexBufferView;
    D3D12_INDEX_BUFFER_VIEW m_indexBufferView;
    uint32_t m_indexCount = 0;
};

class GPUCamera
{
public:
    GPUCamera();
    
    // Camera manipulation methods
    void setPosition(const float3& position);
    void setLookAt(const float3& target);
    void setFOV(float fovInDegrees);
    void setAspectRatio(float aspectRatio);
    
    float3 getPosition() const;
    float3 getDirection() const;
    float getFOV() const;
    
    // Matrix getters
    DirectX::XMMATRIX getViewMatrix() const;
    DirectX::XMMATRIX getProjectionMatrix() const;
    
private:
    float3 m_position = float3(0.0f, 0.0f, -5.0f);
    float3 m_target = float3(0.0f, 0.0f, 0.0f);
    float3 m_up = float3(0.0f, 1.0f, 0.0f);
    float m_fov = 45.0f;
    float m_aspectRatio = 16.0f / 9.0f;
    float m_nearPlane = 0.1f;
    float m_farPlane = 1000.0f;
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

private:
    bool initDirectX();
    
    HWND m_hwnd;
    uint32_t m_width;
    uint32_t m_height;
    std::unique_ptr<UIState> m_uiState;
    
    Microsoft::WRL::ComPtr<ID3D12Device> m_device;
    Microsoft::WRL::ComPtr<IDXGISwapChain4> m_swapChain;
    std::shared_ptr<GPUQueue> m_gpuQueue;
};

class GPUWorld
{
public:
    GPUWorld(std::shared_ptr<Window> pWindow);
    ~GPUWorld();

    // Mesh management
    std::shared_ptr<GPUMesh> createMesh();
    void addMesh(std::shared_ptr<GPUMesh> mesh);
    void removeMesh(std::shared_ptr<GPUMesh> mesh);
    
    // Camera management
    std::shared_ptr<GPUCamera> getCamera();
    void setCamera(std::shared_ptr<GPUCamera> camera);
    
    // Rendering
    void drawMeshesIntoWindow();

private:
    void initializeRenderResources();
    
    struct TransformBuffer {
        DirectX::XMMATRIX World;
        DirectX::XMMATRIX View;
        DirectX::XMMATRIX Projection;
    };
    
    std::shared_ptr<Window> m_pWindow;
    std::vector<std::shared_ptr<GPUMesh>> m_pMeshes;
    std::shared_ptr<GPUCamera> m_pCamera;
    
    // DirectX rendering resources
    Microsoft::WRL::ComPtr<ID3D12RootSignature> m_rootSignature;
    Microsoft::WRL::ComPtr<ID3D12PipelineState> m_pipelineState;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_transformBufferResource;
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_cbvHeap;
    TransformBuffer m_transformBufferMatrix;
    UINT8* m_transformBufferData = nullptr;
};
