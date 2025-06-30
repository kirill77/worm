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