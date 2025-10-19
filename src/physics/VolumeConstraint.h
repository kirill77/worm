#pragma once

#include <vector>
#include <memory>
#include "PhysicsMesh.h"
#include "PhysicsConstraints.h"

// XPBD-style volume constraint operating on PhysicsMesh face view
class VolumeConstraintXPBD : public IConstraint
{
public:
    VolumeConstraintXPBD(PhysicsMesh& body, double targetVolume, double compliance = 0.0)
        : m_body(body)
        , m_targetVolume(targetVolume)
        , m_compliance(compliance)
        , m_lambda(0.0)
    {
    }

    void setTargetVolume(double v) { m_targetVolume = v; }
    double getTargetVolume() const { return m_targetVolume; }
    void setCompliance(double c) { m_compliance = c; }
    double getCompliance() const { return m_compliance; }

    // Project positions to satisfy volume constraint (soft if compliance > 0)
    void project(double dt) override;

    // Utility: compute signed volume using faces and node positions
    double computeSignedVolume() const;

private:
    PhysicsMesh& m_body;
    double m_targetVolume;
    double m_compliance;      // XPBD compliance (0 for hard), in units of 1/stiffness
    double m_lambda;          // XPBD Lagrange multiplier accumulator
};


