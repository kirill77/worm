#include "CameraTransition.h"
#include "visualization/gpu/GPUCamera.h"
#include <algorithm>
#include "geometry/vectors/quat.h"

CameraTransition::CameraTransition(std::shared_ptr<GPUCamera> currentCamera,
                                   std::shared_ptr<GPUCamera> targetCamera,
                                   float durationSec,
                                   const box3& focusBox)
    : m_currentCamera(currentCamera)
    , m_targetCamera(std::move(targetCamera))
    , m_durationSec(std::max(0.0f, durationSec))
    , m_focusBox(focusBox)
{
    if (auto cam = m_currentCamera.lock())
    {
        m_startTransform = cam->getCameraTransform();
        m_startPos = cam->getPosition();
        m_startDir = cam->getDirection();
        m_startUp = cam->getUp();
    }
    if (m_targetCamera)
    {
        m_endTransform = m_targetCamera->getCameraTransform();
        m_endPos = m_targetCamera->getPosition();
        m_endDir = m_targetCamera->getDirection();
        m_endUp = m_targetCamera->getUp();
    }

    if (m_durationSec == 0.0f)
    {
        // Snap immediately
        if (auto cam = m_currentCamera.lock())
        {
            cam->setCameraTransform(m_endTransform);
        }
        m_finished = true;
    }
}

static inline quaternion<float> orientationFromDirUp(const float3& dir, const float3& up)
{
    // Build orthonormal basis
    float3 f = normalize(dir);
    float3 r = normalize(cross(up, f));
    float3 u = normalize(cross(f, r));

    float3x3 m(r, u, f); // rows

    // Build affine from matrix and extract quaternion using helper
    affine3 a(m, float3::zero());
    quaternion<float> q;
    vector<float,3> tmpT, tmpS;
    decomposeAffine(a, &tmpT, &q, &tmpS);
    return normalize(q);
}

bool CameraTransition::update(float deltaTimeSec)
{
    if (m_finished)
        return false;

    m_elapsedSec += std::max(0.0f, deltaTimeSec);
    float t = (m_durationSec > 0.0f) ? (m_elapsedSec / m_durationSec) : 1.0f;
    t = std::max(0.0f, std::min(1.0f, t));

    // Smoothstep for ease-in-out
    float u = smoothstep01(t);

    // Interpolate position linearly
    float3 pos = m_startPos * (1.0f - u) + m_endPos * u;

    // Interpolate orientation using quaternions
    quaternion<float> q0 = orientationFromDirUp(m_startDir, m_startUp);
    quaternion<float> q1 = orientationFromDirUp(m_endDir, m_endUp);
    quaternion<float> q = slerp(q0, q1, u);
    q = normalize(q);
    float3x3 rot = q.toMatrix();
    float3 right = rot.row0;
    float3 up = rot.row1;
    float3 dir = rot.row2;

    // Build transform
    affine3 xform;
    xform.m_linear = float3x3(right, up, dir);
    xform.m_translation = pos;

    if (auto cam = m_currentCamera.lock())
    {
        cam->setCameraTransform(xform);
    }

    if (t >= 1.0f)
    {
        m_finished = true;
        return false;
    }
    return true;
}


