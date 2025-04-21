#pragma once

#include "../math/vector.h"
#include <DirectXMath.h>

class GPUCamera
{
public:
    GPUCamera();
    
    // Camera manipulation methods
    void setPosition(const float3& pos);
    void setDirection(const float3& dir);
    void setUp(const float3& up);
    void setFOV(float fov);
    void setAspectRatio(float aspectRatio);
    
    float3 getPosition() const;
    float3 getDirection() const;
    float3 getUp() const { return m_up; }
    float3 getRight() const;
    float getFOV() const;
    
    // Matrix getters
    DirectX::XMMATRIX getViewMatrix() const;
    DirectX::XMMATRIX getProjectionMatrix() const;
    
private:
    float3 m_position = float3(0.0f, 0.0f, -5.0f);
    float3 m_direction = float3(0.0f, 0.0f, 1.0f);  // Default looking along positive Z
    float3 m_up = float3(0.0f, 1.0f, 0.0f);
    float m_fov = 45.0f;
    float m_aspectRatio = 16.0f / 9.0f;
    float m_nearPlane = 0.1f;
    float m_farPlane = 1000.0f;
}; 