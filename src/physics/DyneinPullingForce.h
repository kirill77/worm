#pragma once

#include "ForceGenerator.h"
#include <vector>
#include <memory>

class PhysCentrosome;

// Force generator for dynein motor pulling on cortex-bound microtubules
// Applies pulling forces to the cortex mesh at microtubule attachment points
class DyneinPullingForce : public IForceGenerator
{
public:
    DyneinPullingForce(
        PhysicsMesh& cortexBody,
        const std::vector<std::shared_ptr<PhysCentrosome>>& centrosomes,
        double pullingForcePerMT);

    void apply(double dt) override;

private:
    PhysicsMesh& m_cortexBody;
    const std::vector<std::shared_ptr<PhysCentrosome>>& m_centrosomes;
    double m_pullingForcePerMT; // Force magnitude per bound microtubule (in newtons or simulation units)
};

