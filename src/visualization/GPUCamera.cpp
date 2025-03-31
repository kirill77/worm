#include "pch.h"
#include "Window.h"
#include "DirectXHelpers.h"

GPUCamera::GPUCamera()
{
    // Default initialization is done in the class definition
}

void GPUCamera::setPosition(const float3& position)
{
    m_position = position;
}

void GPUCamera::setLookAt(const float3& target)
{
    m_target = target;
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
    // Calculate direction vector from position to target
    return normalize(m_target - m_position);
}

float GPUCamera::getFOV() const
{
    return m_fov;
}

DirectX::XMMATRIX GPUCamera::getViewMatrix() const
{
    // Convert float3 to XMVECTOR
    DirectX::XMVECTOR eyePosition = DirectX::XMVectorSet(m_position.x, m_position.y, m_position.z, 1.0f);
    DirectX::XMVECTOR focusPosition = DirectX::XMVectorSet(m_target.x, m_target.y, m_target.z, 1.0f);
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