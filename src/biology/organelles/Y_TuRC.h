#pragma once

#include <memory>
#include <vector>
#include "geometry/vectors/vector.h"
#include "chemistry/molecules/Molecule.h"
#include "physics/PhysMicrotubule.h"
#include "geometry/mesh/MeshLocation.h"

struct Centrosome;
#include "Cortex.h" // needed for nested Cortex::IntersectionResult type
class Medium;

// represents Gamma-TuRC (gamma-tubulin ring complex) - the place where a microtubule nucleates
struct Y_TuRC : public PhysMicrotubule
{
    Y_TuRC(std::weak_ptr<Centrosome> pCentrosome);

    // Update microtubule lifecycle and dynamics (simple dynamic instability)
    // centrosomeCellPos: cell-space position of the centrosome center (µm, cortex-centered)
    // pCortex: cortex organelle for geometry queries
    void update(double dtSec, const float3& centrosomeCellPos, const std::shared_ptr<Cortex>& pCortex, Medium& internalMedium);
    
    // Y_TuRC-specific accessor
    float3 getTipDirection() const { return m_vTipDir; }

private:
    std::weak_ptr<Centrosome> m_pCentrosome;
    Species m_species = Species::GENERIC;

    float3 m_vTipDir; // current tip direction (normalized)

    // Microtubule refractory and contact state
    float m_mtRefractorySec = 0.0f; // wait time before re-nucleation after disassembly
    bool  m_mtContactCortex = false;   // whether tip is in contact with cortex
    
    // Cortical binding state
    double m_bindingStrength = 0.0;   // dynein-dependent binding strength for unbinding kinetics
    double m_bindingTime = 0.0;       // how long MT has been bound to cortex

    // Helper methods
    float getCapLengthMicroM() const;
    
    // Cortical binding state management 
    void bindToCortex(const MeshLocation& location, double dyneinConc);
    void unbindFromCortex();
    bool shouldUnbind(double dtSec, const std::function<float()>& rand01) const;
    
    // Enhanced binding with cached intersection details
    void attemptCorticalBindingWithIntersection(const Cortex::CortexRay& intersection,
        const std::shared_ptr<Cortex>& pCortex, Medium& internalMedium, 
        double dtSec, const std::function<float()>& rand01);
};

