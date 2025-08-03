#include "pch.h"
#include "Window.h"
#include "DirectXHelpers.h"
#include "GPUCamera.h"

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
    // The linear part contains the camera's orientation as column vectors
    affine3 transform;
    transform.m_linear = float3x3::from_cols(right, up, forward);
    transform.m_translation = m_position;
    
    return transform;
}

void GPUCamera::setCameraTransform(const affine3& transform)
{
    // Extract position from translation
    m_position = transform.m_translation;
    
    // Extract orientation vectors from linear transform (column vectors)
    float3 right = transform.m_linear.col(0);
    float3 up = transform.m_linear.col(1);
    float3 forward = transform.m_linear.col(2);
    
    // Normalize to ensure orthonormal basis
    m_direction = normalize(forward);
    m_up = normalize(up);
    
    // Ensure right vector is consistent (cross product should give orthogonal result)
    // We trust the up and direction vectors more than the right vector
    float3 computedRight = normalize(cross(m_up, m_direction));
    
    // Re-orthogonalize: if the input wasn't perfectly orthonormal, fix it
    m_up = normalize(cross(m_direction, computedRight));
} 