#pragma once

#include <memory>
#include "geometry/vectors/vector.h"

struct Centrosome;
class Cortex;

// represents Gamma-TuRC (gamma-tubulin ring complex) - the place where a microtubule nucleates
struct Y_TuRC
{
    Y_TuRC(std::weak_ptr<Centrosome> pCentrosome);

    const float3& getDirection() { return m_vDir; }
    const float3& getPosition() { return m_vPosMicroM; }

    // Update microtubule lifecycle and dynamics (simple dynamic instability)
    // centrosomeWorldPos: world-space position of the centrosome center (µm)
    // pCortex: cortex organelle for geometry queries
    void update(double dtSec, const float3& centrosomeWorldPos, const std::shared_ptr<Cortex>& pCortex);
    // Accessors for MT visualization (optional)
    float getMTLengthMicroM() const { return m_mtLengthMicroM; }
    bool  hasActiveMT() const { return m_mtLengthMicroM > 0.0f; }

private:
    std::weak_ptr<Centrosome> m_pCentrosome;

    float3 m_vDir; // normalized direction
    float3 m_vPosMicroM; // position in micro-meters in respect to centrosome center

    // Microtubule dynamic state (inline, without a separate class)
    enum class MTState { Growing, Shrinking };
    MTState m_mtState = MTState::Growing;
    float m_mtLengthMicroM = 0.0f; // current length in µm
    float m_mtRefractorySec = 0.0f; // wait time before re-nucleation after disassembly

    // microtubule grows by attaching alpha and beta tubulins
    uint32_t m_nAlphaTubulins = 0;
    uint32_t m_nBetaTubulins = 0;
};

