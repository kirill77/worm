#pragma once

#include <vector>
#include <memory>
#include "geometry/vectors/vector.h"
#include "physics/PhysicsMesh.h"

/**
 * Interface for force generators acting on a mesh-based soft body.
 * Each force is bound to specific bodies at construction time.
 */
class IForceGenerator
{
public:
    virtual ~IForceGenerator() = default;

    // Apply forces to the associated body
    virtual void apply() = 0;
};

/** Edge-aligned Hookean springs for each mesh edge */
class EdgeSpringForce : public IForceGenerator
{
public:
    EdgeSpringForce(PhysicsMesh& body, double springConstant);

    void apply() override;

private:
    PhysicsMesh& m_body;
    double m_springConstant;
    std::vector<double> m_edgeRestLengths;
};

/** Edge-aligned relative-velocity damping for each mesh edge */
class EdgeDampingForce : public IForceGenerator
{
public:
    EdgeDampingForce(PhysicsMesh& body, double dampingCoeff)
        : m_body(body)
        , m_dampingCoeff(dampingCoeff) {
    }

    void apply() override;

private:
    PhysicsMesh& m_body;
    double m_dampingCoeff;
};


