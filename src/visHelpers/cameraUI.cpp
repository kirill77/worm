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

void CameraUI::moveLeft(float fDistance)
{
    if (!m_pCamera)
        return;

    // Get current camera position and direction
    float3 cameraPos = m_pCamera->getPosition();
    float3 direction = m_pCamera->getDirection();
    direction = normalize(direction);
    
    // Calculate left vector using camera's up vector
    float3 vRight = m_pCamera->getRight();
    
    // Move camera left by the given distance
    float3 newPosition = cameraPos - vRight * fDistance;
    
    // Update camera position
    m_pCamera->setPosition(newPosition);
}

void CameraUI::notifyNewUIState(const UIState& uiState)
{
    if (!m_pCamera)
        return;

    // Handle Ctrl+A to fit world box in view
    if (uiState.isButtonOrKeyPressed(VK_CONTROL) &&
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
    if (uiState.isButtonOrKeyPressed(VK_LBUTTON)) // Left mouse button for world rotation
    {
        // Calculate rotation based on mouse movement
        float2 currentMousePos = uiState.getMousePosition();
        float2 prevMousePos = m_prevUIState.getMousePosition();
        float deltaX = currentMousePos.x - prevMousePos.x;
        float deltaY = currentMousePos.y - prevMousePos.y;
        
        // Get current camera position and direction
        float3 cameraPos = m_pCamera->getPosition();
        float3 direction = m_pCamera->getDirection();
        float3 up = m_pCamera->getUp();
        float3 right = m_pCamera->getRight();
        float distance = length(direction);
        direction = normalize(direction);
        
        // Calculate the world center (use origin if world box is empty)
        float3 worldCenter = m_worldBox.isempty() ? float3(0.0f, 0.0f, 0.0f) : m_worldBox.center();
        
        // Calculate the vector from camera to world center
        float3 centerToCamera = worldCenter - cameraPos;
        
        // Calculate rotation angles
        float yawAngle = deltaX * m_rotationSpeed * 0.01f;
        float pitchAngle = deltaY * m_rotationSpeed * 0.01f;
        
        // Create rotation quaternions
        DirectX::XMVECTOR yawQuat = DirectX::XMQuaternionRotationAxis(
            DirectX::XMVectorSet(up.x, up.y, up.z, 0.0f),
            yawAngle
        );
        DirectX::XMVECTOR pitchQuat = DirectX::XMQuaternionRotationAxis(
            DirectX::XMVectorSet(right.x, right.y, right.z, 0.0f),
            pitchAngle
        );
        
        // Combine rotations
        DirectX::XMVECTOR rotationQuat = DirectX::XMQuaternionMultiply(pitchQuat, yawQuat);
        
        // Rotate the center-to-camera vector
        DirectX::XMVECTOR centerToCameraVec = DirectX::XMVectorSet(centerToCamera.x, centerToCamera.y, centerToCamera.z, 0.0f);
        DirectX::XMVECTOR rotatedVector = DirectX::XMVector3Rotate(centerToCameraVec, rotationQuat);
        
        // Calculate new camera position and direction
        DirectX::XMFLOAT3 newCenterToCamera;
        DirectX::XMStoreFloat3(&newCenterToCamera, rotatedVector);
        float3 newCameraPos = worldCenter - float3(newCenterToCamera.x, newCenterToCamera.y, newCenterToCamera.z);
        float3 newDirection = normalize(worldCenter - newCameraPos);
        
        // Rotate the up vector
        DirectX::XMVECTOR upVec = DirectX::XMVectorSet(up.x, up.y, up.z, 0.0f);
        DirectX::XMVECTOR rotatedUp = DirectX::XMVector3Rotate(upVec, rotationQuat);
        DirectX::XMFLOAT3 newUp;
        DirectX::XMStoreFloat3(&newUp, rotatedUp);
        
        // Update camera position, direction, and up vector
        m_pCamera->setPosition(newCameraPos);
        m_pCamera->setDirection(newDirection);
        m_pCamera->setUp(float3(newUp.x, newUp.y, newUp.z));
    }

    // rotation of camera around camera center
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
        float3 right = m_pCamera->getRight();
        float3 up = m_pCamera->getUp();
        
        // Apply rotation around up vector (yaw)
        float yawAngle = -deltaX * m_rotationSpeed * 0.01f;
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

    if (uiState.isButtonOrKeyPressed('W'))
    {
        moveForward(calculateMoveSpeed());
    }
    if (uiState.isButtonOrKeyPressed('S'))
    {
        moveForward(-calculateMoveSpeed());
    }
    if (uiState.isButtonOrKeyPressed('A'))
    {
         moveLeft(calculateMoveSpeed());
    }
    if (uiState.isButtonOrKeyPressed('D'))
    {
        moveLeft(-calculateMoveSpeed());
    }
    
    // Store the current UI state for the next frame
    m_prevUIState = uiState;
}
