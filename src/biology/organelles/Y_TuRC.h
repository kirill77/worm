#pragma once

#include <memory>
#include <vector>
#include "geometry/vectors/vector.h"
#include "chemistry/molecules/Molecule.h"

struct Centrosome;
#include "Cortex.h" // needed for nested Cortex::IntersectionResult type
class Medium;

// represents Gamma-TuRC (gamma-tubulin ring complex) - the place where a microtubule nucleates
struct Y_TuRC : public std::enable_shared_from_this<Y_TuRC>
{
    // Microtubule dynamic state (inline, without a separate class)
    enum class MTState { Growing, Shrinking, Bound };

    Y_TuRC(std::weak_ptr<Centrosome> pCentrosome);

    const float3& getOrigin() const { return m_mtSegmentPoints[0]; }
    const std::vector<float3>& getSegmentPoints() const { return m_mtSegmentPoints; }

    // Update microtubule lifecycle and dynamics (simple dynamic instability)
    // centrosomeWorldPos: world-space position of the centrosome center (µm)
    // pCortex: cortex organelle for geometry queries
    void update(double dtSec, const float3& centrosomeWorldPos, const std::shared_ptr<Cortex>& pCortex, Medium& internalMedium);
    // Accessors for MT visualization (optional)
    float getMTLengthMicroM() const;
    bool  hasActiveMT() const { return m_mtSegmentPoints.size() >= 2; }
    float3 getTipPosition() const { return m_mtSegmentPoints.back(); }
    float3 getTipDirection() const { return m_vTipDir; }
    
    // Cortical binding accessors
    bool isBound() const { return m_mtState == MTState::Bound; }
    MTState getState() const { return m_mtState; }

private:
    std::weak_ptr<Centrosome> m_pCentrosome;
    Species m_species = Species::GENERIC;

    float3 m_vTipDir; // current tip direction (normalized)
    std::vector<float3> m_mtSegmentPoints; // array of points: [0] = nucleation origin, [1..n] = segment endpoints (always size >= 2)

    // Microtubule dynamic state
    MTState m_mtState = MTState::Growing;
    float m_mtRefractorySec = 0.0f; // wait time before re-nucleation after disassembly
    bool  m_mtContactCortex = false;   // whether tip is in contact with cortex
    
    // Cortical binding state
    double m_bindingStrength = 0.0;   // dynein-dependent binding strength for unbinding kinetics
    double m_bindingTime = 0.0;       // how long MT has been bound to cortex

    // Helper methods
    float getLastSegmentLength() const;
    float getCapLengthMicroM() const;
    
    // Cortical binding state management 
    void bindToCortex(double dyneinConc);
    void unbindFromCortex();
    bool shouldUnbind(double dtSec, const std::function<float()>& rand01) const;
    
    // Enhanced binding with cached intersection details
    void attemptCorticalBindingWithIntersection(const Cortex::CortexRay& intersection,
        const std::shared_ptr<Cortex>& pCortex, Medium& internalMedium, 
        double dtSec, const std::function<float()>& rand01);
};

