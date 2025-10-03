#pragma once

#include <memory>
#include "geometry/vectors/vector.h"
#include "chemistry/molecules/Molecule.h"

struct Centrosome;
class Cortex;
class Medium;

// represents Gamma-TuRC (gamma-tubulin ring complex) - the place where a microtubule nucleates
struct Y_TuRC
{
    Y_TuRC(std::weak_ptr<Centrosome> pCentrosome);

    const float3& getDirection() { return m_vDir; }
    const float3& getPosition() { return m_vPosMicroM; }

    // Update microtubule lifecycle and dynamics (simple dynamic instability)
    // centrosomeWorldPos: world-space position of the centrosome center (µm)
    // pCortex: cortex organelle for geometry queries
    void update(double dtSec, const float3& centrosomeWorldPos, const std::shared_ptr<Cortex>& pCortex, Medium& internalMedium);
    // Accessors for MT visualization (optional)
    float getMTLengthMicroM() const { return m_mtLengthMicroM; }
    bool  hasActiveMT() const { return m_mtLengthMicroM > 0.0f; }

private:
    std::weak_ptr<Centrosome> m_pCentrosome;
    Species m_species = Species::GENERIC;

    float3 m_vDir; // normalized direction
    float3 m_vPosMicroM; // position in micro-meters in respect to centrosome center

    // Microtubule dynamic state (inline, without a separate class)
    enum class MTState { Growing, Shrinking };
    MTState m_mtState = MTState::Growing;
    float m_mtLengthMicroM = 0.0f; // current length in µm
    float m_mtRefractorySec = 0.0f; // wait time before re-nucleation after disassembly
    float m_mtCapLengthMicroM = 0.0f; // GTP-cap proxy length in µm
    bool  m_mtContactCortex = false;   // whether tip is in contact with cortex

    // microtubule grows by attaching alpha and beta tubulins
    uint32_t m_nAlphaTubulins = 0;
    uint32_t m_nBetaTubulins = 0;
};

