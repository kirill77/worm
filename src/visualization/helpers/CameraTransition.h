#pragma once

#include <memory>
#include "geometry/vectors/affine.h"
#include "geometry/vectors/box.h"

class GPUCamera;

// Smoothly transitions a camera from its current transform to a target camera's transform
// over a specified duration. Call update(dtSec) every frame until isFinished() is true.
class CameraTransition
{
public:
    CameraTransition(std::shared_ptr<GPUCamera> currentCamera,
                     std::shared_ptr<GPUCamera> targetCamera,
                     float durationSec,
                     const box3& focusBox);

    // Advance the transition. Returns true while transition is in progress, false when finished.
    bool update(float deltaTimeSec);

    bool isFinished() const { return m_finished; }
    void cancel() { m_finished = true; }
    
    // Get the focus box that this transition is targeting
    const box3& getFocusBox() const { return m_focusBox; }

private:
    std::weak_ptr<GPUCamera> m_currentCamera;
    std::shared_ptr<GPUCamera> m_targetCamera;

    // Cached endpoints
    affine3 m_startTransform;
    affine3 m_endTransform;

    // Decomposed endpoints (for stable interpolation)
    float3 m_startPos;
    float3 m_endPos;
    float3 m_startDir;
    float3 m_endDir;
    float3 m_startUp;
    float3 m_endUp;

    float m_durationSec = 0.0f;
    float m_elapsedSec = 0.0f;
    bool m_finished = false;
    
    // The bounding box that the camera is focusing on
    box3 m_focusBox;

    static inline float smoothstep01(float t)
    {
        t = std::max(0.0f, std::min(1.0f, t));
        return t * t * (3.0f - 2.0f * t);
    }
};


