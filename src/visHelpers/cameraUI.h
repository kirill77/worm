#pragma once

#include "visualization/Window.h"
#include <memory>

class GPUCamera;

class CameraUI
{
public:
    CameraUI();
    ~CameraUI();

    // Attach a camera to control
    void attachToCamera(std::shared_ptr<GPUCamera> pCamera);
    
    // Update camera based on UI state changes
    void notifyNewUIState(const UIState& uiState);

private:
    // Camera movement parameters
    float m_rotationSpeed = 0.01f;
    float m_panSpeed = 0.01f;
    float m_zoomSpeed = 0.1f;
    
    // Camera reference
    std::shared_ptr<GPUCamera> m_pCamera;
    
    // Previous UI state for delta calculations
    UIState m_prevUIState;
    
    // Helper methods
    void handleRotation(const float2& mouseDelta);
    void handlePanning(const float2& mouseDelta);
    void handleZooming(float scrollDelta);
};

