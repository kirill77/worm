#include "cameraUI.h"
#include "visualization/GPUCamera.h"
#include <cmath>

CameraUI::CameraUI()
    : m_rotationSpeed(0.5f)
    , m_panSpeed(0.01f)
    , m_zoomSpeed(0.1f)
    , m_pCamera(nullptr)
{
}

CameraUI::~CameraUI()
{
}

void CameraUI::attachToCamera(std::shared_ptr<GPUCamera> camera)
{
    m_pCamera = camera;
}

float CameraUI::calculateMoveSpeed() const
{
    float moveSpeed = 0.1f;
    if (!m_worldBox.isempty())
    {
        float3 boxDiagonal = m_worldBox.diagonal();
        float maxDimension = std::max(std::max(boxDiagonal.x, boxDiagonal.y), boxDiagonal.z);
        moveSpeed = maxDimension * 0.1f; // Move 10% of the world size per frame
    }
    return moveSpeed;
}

void CameraUI::moveForward(float fDistance)
{
    if (!m_pCamera)
        return;

    // Get current camera position and direction
    float3 cameraPos = m_pCamera->getPosition();
    float3 direction = m_pCamera->getDirection();
    direction = normalize(direction);
    
    // Move camera forward by the given distance
    float3 newPosition = cameraPos + direction * fDistance;
    
    // Update camera position
    m_pCamera->setPosition(newPosition);
}

void CameraUI::notifyNewUIState(const UIState& uiState)
{
    if (!m_pCamera)
        return;

    // Handle Ctrl+A to fit world box in view
    if (uiState.getButtonOrKeyPressCount(VK_CONTROL) &&
        uiState.getButtonOrKeyPressCount('A') > m_prevUIState.getButtonOrKeyPressCount('A'))
    {
        if (!m_worldBox.isempty())
        {
            // Set FOV to 30 degrees
            m_pCamera->setFOV(30.0f);
            
            // Calculate the center of the world box
            float3 boxCenter = m_worldBox.center();
            
            // Calculate the diagonal of the world box
            float3 boxDiagonal = m_worldBox.diagonal();
            
            // Calculate the maximum dimension of the box
            float maxDimension = std::max(std::max(boxDiagonal.x, boxDiagonal.y), boxDiagonal.z);
            
            // Calculate the distance needed to fit the box in view
            // For a 30-degree FOV, we need to be at a distance of maxDimension / (2 * tan(15 degrees))
            float fovRadians = 30.0f * 3.14159265359f / 180.0f;
            float distance = maxDimension / (2.0f * std::tan(fovRadians / 2.0f));
            
            // Add a small margin to ensure the box is fully visible
            distance *= 1.1f;
            
            // Calculate the new camera position
            // Position the camera at the center of the box, offset by the calculated distance
            // along the negative z-axis (assuming the camera looks along the positive z-axis)
            float3 newPosition = boxCenter;
            newPosition.z -= distance;
            
            // Set the camera position and target
            m_pCamera->setPosition(newPosition);
            m_pCamera->setDirection(boxCenter - newPosition);
            
            return; // Skip other camera controls when fitting to view
        }
    }

    // Handle rotation
    if (uiState.isButtonOrKeyPressed(VK_RBUTTON)) // Right mouse button
    {
        // Calculate rotation based on mouse movement
        float2 currentMousePos = uiState.getMousePosition();
        float2 prevMousePos = m_prevUIState.getMousePosition();
        float deltaX = currentMousePos.x - prevMousePos.x;
        float deltaY = currentMousePos.y - prevMousePos.y;
        
        // Get current camera position and direction
        float3 cameraPos = m_pCamera->getPosition();
        float3 direction = m_pCamera->getDirection();
        float distance = length(direction);
        direction = normalize(direction);
        
        // Calculate right and up vectors
        float3 right = normalize(cross(direction, float3(0.0f, 1.0f, 0.0f)));
        float3 up = normalize(cross(right, direction));
        
        // Apply rotation around up vector (yaw)
        float yawAngle = deltaX * m_rotationSpeed * 0.01f;
        float cosYaw = std::cos(yawAngle);
        float sinYaw = std::sin(yawAngle);
        
        float3 newRight = right * cosYaw + direction * sinYaw;
        float3 newDirection = -right * sinYaw + direction * cosYaw;
        
        // Apply rotation around right vector (pitch)
        float pitchAngle = deltaY * m_rotationSpeed * 0.01f;
        float cosPitch = std::cos(pitchAngle);
        float sinPitch = std::sin(pitchAngle);
        
        float3 finalDirection = newDirection * cosPitch - up * sinPitch;
        
        // Update camera target
        m_pCamera->setDirection(finalDirection);
    }
    
    // Handle panning
    if (uiState.isButtonOrKeyPressed(VK_LBUTTON))
    {
        // Calculate pan offset based on mouse movement
        float2 currentMousePos = uiState.getMousePosition();
        float2 prevMousePos = m_prevUIState.getMousePosition();
        float deltaX = currentMousePos.x - prevMousePos.x;
        float deltaY = currentMousePos.y - prevMousePos.y;
        
        // Get current camera position and direction
        float3 cameraPos = m_pCamera->getPosition();
        float3 direction = m_pCamera->getDirection();
        float distance = length(direction);
        direction = normalize(direction);
        
        // Calculate right and up vectors
        float3 right = normalize(cross(direction, float3(0.0f, 1.0f, 0.0f)));
        float3 up = normalize(cross(right, direction));
        
        // Calculate the pan offset
        float3 panOffset = -right * deltaX * m_panSpeed * distance - up * deltaY * m_panSpeed * distance;
        
        // Update camera position and target
        m_pCamera->setPosition(cameraPos + panOffset);
    }
    
    // Handle zooming with mouse wheel
    float scrollDelta = uiState.getScrollWheelState();
    if (std::abs(scrollDelta) > 0.0f)
    {
        // Get current camera position and direction
        float3 cameraPos = m_pCamera->getPosition();
        float3 direction = m_pCamera->getDirection();
        float distance = length(direction);
        direction = normalize(direction);
        
        // Calculate the zoom factor
        float zoomFactor = 1.0f + scrollDelta * m_zoomSpeed * 0.01f;
        
        // Calculate the new distance
        float newDistance = distance * zoomFactor;
        
        // Ensure the new distance is not too small
        newDistance = std::max(newDistance, 0.1f);
        
        // Calculate the new camera position
        float3 newPosition = cameraPos + direction * (newDistance - distance);
        
        // Update camera position
        m_pCamera->setPosition(newPosition);
    }

    // Handle forward movement with 'w' key
    if (uiState.isButtonOrKeyPressed('W'))
    {
        // Move camera forward
        moveForward(calculateMoveSpeed());
    }
    // Handle backward movement with 's' key
    if (uiState.isButtonOrKeyPressed('S'))
    {
        // Move camera forward
        moveForward(-calculateMoveSpeed());
    }
    
    // Store the current UI state for the next frame
    m_prevUIState = uiState;
}
