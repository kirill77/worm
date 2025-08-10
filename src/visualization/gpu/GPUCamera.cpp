#include "pch.h"
#include "Window.h"
#include "DirectXHelpers.h"
#include "GPUCamera.h"
#include <algorithm>
#include <cmath>

GPUCamera::GPUCamera()
{
    // Default initialization is done in the class definition
}

void GPUCamera::setPosition(const float3& position)
{
    m_position = position;
}

void GPUCamera::setDirection(const float3& direction)
{
    m_direction = normalize(direction);
}

void GPUCamera::setUp(const float3& up)
{
    m_up = normalize(up);
}

void GPUCamera::setFOV(float fovInDegrees)
{
    m_fov = fovInDegrees;
}

void GPUCamera::setAspectRatio(float aspectRatio)
{
    m_aspectRatio = aspectRatio;
}

float3 GPUCamera::getPosition() const
{
    return m_position;
}

float3 GPUCamera::getDirection() const
{
    return m_direction;
}

float3 GPUCamera::getRight() const
{
    return normalize(cross(m_up, m_direction));
}

float GPUCamera::getFOV() const
{
    return m_fov;
}

DirectX::XMMATRIX GPUCamera::getViewMatrix() const
{
    // Convert float3 to XMVECTOR
    DirectX::XMVECTOR eyePosition = DirectX::XMVectorSet(m_position.x, m_position.y, m_position.z, 1.0f);
    DirectX::XMVECTOR focusPosition = DirectX::XMVectorSet(
        m_position.x + m_direction.x,
        m_position.y + m_direction.y,
        m_position.z + m_direction.z,
        1.0f);
    DirectX::XMVECTOR upDirection = DirectX::XMVectorSet(m_up.x, m_up.y, m_up.z, 0.0f);
    
    // Create view matrix
    return DirectX::XMMatrixLookAtLH(eyePosition, focusPosition, upDirection);
}

DirectX::XMMATRIX GPUCamera::getProjectionMatrix() const
{
    // Create projection matrix with specified FOV
    float fovRadians = DirectX::XMConvertToRadians(m_fov);
    return DirectX::XMMatrixPerspectiveFovLH(fovRadians, m_aspectRatio, m_nearPlane, m_farPlane);
}

affine3 GPUCamera::getCameraTransform() const
{
    // Calculate camera basis vectors
    float3 right = getRight();
    float3 up = m_up;
    float3 forward = m_direction;
    
    // Create camera-to-world transform
    // The linear part contains the camera's orientation as row vectors
    affine3 transform;
    transform.m_linear = float3x3(right, up, forward);
    transform.m_translation = m_position;
    
    return transform;
}

bool GPUCamera::fitBoxToView(const box3& boxToFit)
{
    if (boxToFit.isempty())
        return false;

    // Set FOV to 30 degrees for a consistent framing
    setFOV(30.0f);

    // Calculate the center and size of the box
    float3 boxCenter = boxToFit.center();
    float3 boxDiagonal = boxToFit.diagonal();
    float maxDimension = std::max({boxDiagonal.x, boxDiagonal.y, boxDiagonal.z});

    // Distance to fit the box given FOV
    float fovRadians = DirectX::XMConvertToRadians(getFOV());
    float halfFov = fovRadians * 0.5f;
    float distance = maxDimension / (2.0f * std::tan(halfFov));

    // Margin to ensure full visibility
    distance *= 1.1f;

    // Position the camera along its current forward direction
    float3 forward = normalize(getDirection());
    if (!std::isfinite(forward.x) || !std::isfinite(forward.y) || !std::isfinite(forward.z))
        forward = float3(0.0f, 0.0f, 1.0f);

    float3 newPosition = boxCenter - forward * distance;

    setPosition(newPosition);
    setDirection(boxCenter - newPosition);
    return true;
}

void GPUCamera::setCameraTransform(const affine3& transform)
{
    // Extract position from translation
    m_position = transform.m_translation;
    
    // Extract orientation vectors from linear transform (row vectors)
    float3 right = transform.m_linear.row0;
    float3 up = transform.m_linear.row1;
    float3 forward = transform.m_linear.row2;
    
    // Normalize to ensure orthonormal basis
    m_direction = normalize(forward);
    m_up = normalize(up);
    
    // Ensure right vector is consistent (cross product should give orthogonal result)
    // We trust the up and direction vectors more than the right vector
    float3 computedRight = normalize(cross(m_up, m_direction));
    
    // Re-orthogonalize: if the input wasn't perfectly orthonormal, fix it
    m_up = normalize(cross(m_direction, computedRight));
} 