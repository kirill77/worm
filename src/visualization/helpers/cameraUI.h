#pragma once

#include <memory>
#include "geometry/vectors/box.h"
#include "visualization/gpu/Window.h"

class GPUCamera;

class CameraUI
{
public:
    CameraUI();
    ~CameraUI();

    void setWorldBox(const box3& worldBox) { m_worldBox = worldBox; }

    // Attach a camera to control
    void attachToCamera(std::shared_ptr<GPUCamera> pCamera);
    
    // Update camera based on UI state changes
    void notifyNewUIState(const UIState& uiState);
    
    // Fit the world box to the camera view
    bool fitWorldBoxToView();

private:
    // Camera movement parameters
    float m_rotationSpeed = 0.01f;
    float m_panSpeed = 0.01f;
    float m_zoomSpeed = 0.1f;

    box3 m_worldBox;
    
    // Camera reference
    std::shared_ptr<GPUCamera> m_pCamera;
    
    // Previous UI state for delta calculations
    UIState m_prevUIState;
    
    // Helper methods
    void handleRotation(const float2& mouseDelta);
    void handlePanning(const float2& mouseDelta);
    void handleZooming(float scrollDelta);

    void moveForward(float fDistance);
    void moveLeft(float fDistance);
    float calculateMoveSpeed() const;
};

