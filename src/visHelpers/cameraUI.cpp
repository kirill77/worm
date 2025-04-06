#include "cameraUI.h"
#include "visualization/GPUCamera.h"
#include "../math/vector.h"
#include <cmath>

CameraUI::CameraUI()
    : m_pCamera(nullptr)
{
}

CameraUI::~CameraUI()
{
}

void CameraUI::attachToCamera(std::shared_ptr<GPUCamera> pCamera)
{
    m_pCamera = pCamera;
}

void CameraUI::notifyNewUIState(const UIState& uiState)
{
    if (!m_pCamera)
        return;

    // Calculate mouse movement delta
    float2 currentMousePos = uiState.getMousePosition();
    float2 mouseDelta = currentMousePos - m_prevUIState.getMousePosition();
    
    // Handle mouse buttons for camera control
    // Left button (VK_LBUTTON = 0x01) for rotation
    if (uiState.getButtonOrKeyPressCount(0x01) > 0)
    {
        handleRotation(mouseDelta);
    }
    
    // Right button (VK_RBUTTON = 0x02) for panning
    if (uiState.getButtonOrKeyPressCount(0x02) > 0)
    {
        handlePanning(mouseDelta);
    }
    
    // Middle button or scroll wheel for zooming
    float scrollDelta = uiState.getScrollWheelState();
    if (std::abs(scrollDelta) > 0.0f)
    {
        handleZooming(scrollDelta);
    }
    
    // Store current UI state for next frame
    m_prevUIState = uiState;
}

void CameraUI::handleRotation(const float2& mouseDelta)
{
    // Get current camera position and direction
    float3 position = m_pCamera->getPosition();
    float3 direction = m_pCamera->getDirection();
    float distance = length(direction);
    direction = normalize(direction);
    
    // Calculate right vector (cross product of direction and up)
    float3 right = normalize(cross(direction, float3(0.0f, 1.0f, 0.0f)));
    
    // Calculate up vector (cross product of right and direction)
    float3 up = normalize(cross(right, direction));
    
    // Apply rotation around up vector (yaw)
    float yawAngle = mouseDelta.x * m_rotationSpeed;
    float cosYaw = std::cos(yawAngle);
    float sinYaw = std::sin(yawAngle);
    
    float3 newRight = right * cosYaw + direction * sinYaw;
    float3 newDirection = -right * sinYaw + direction * cosYaw;
    
    // Apply rotation around right vector (pitch)
    float pitchAngle = mouseDelta.y * m_rotationSpeed;
    float cosPitch = std::cos(pitchAngle);
    float sinPitch = std::sin(pitchAngle);
    
    float3 finalDirection = newDirection * cosPitch - up * sinPitch;
    
    // Update camera target
    float3 newTarget = position + finalDirection * distance;
    m_pCamera->setLookAt(newTarget);
}

void CameraUI::handlePanning(const float2& mouseDelta)
{
    // Get current camera position and direction
    float3 position = m_pCamera->getPosition();
    float3 direction = m_pCamera->getDirection();
    float distance = length(direction);
    direction = normalize(direction);
    
    // Calculate right and up vectors
    float3 right = normalize(cross(direction, float3(0.0f, 1.0f, 0.0f)));
    float3 up = normalize(cross(right, direction));
    
    // Calculate pan offset
    float3 panOffset = -right * mouseDelta.x * m_panSpeed * distance + 
                       up * mouseDelta.y * m_panSpeed * distance;
    
    // Update camera position and target
    m_pCamera->setPosition(position + panOffset);
    float3 newTarget = position + direction * distance + panOffset;
    m_pCamera->setLookAt(newTarget);
}

void CameraUI::handleZooming(float scrollDelta)
{
    // Get current camera position and direction
    float3 position = m_pCamera->getPosition();
    float3 direction = m_pCamera->getDirection();
    float distance = length(direction);
    direction = normalize(direction);
    
    // Calculate new distance with zoom speed
    float newDistance = distance * (1.0f - scrollDelta * m_zoomSpeed);
    
    // Clamp distance to prevent getting too close or too far
    newDistance = std::max(0.1f, std::min(newDistance, 1000.0f));
    
    // Update camera position
    float3 newPosition = position + direction * (newDistance - distance);
    m_pCamera->setPosition(newPosition);
}
