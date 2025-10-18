#pragma once

#include <vector>
#include <memory>
#include "geometry/vectors/vector.h"
#include "physics/BodyInterfaces.h"

/**
 * Interface for force generators acting on a mesh-based soft body.
 * Implementations write forces into the provided outForces buffer.
 */
class IForceGenerator
{
public:
    virtual ~IForceGenerator() = default;

    // Apply forces to a generic body via node/edge views
    virtual void apply(IFaceBody& body, double dt) = 0;
};

/** Edge-aligned Hookean springs for each mesh edge */
class EdgeSpringForce : public IForceGenerator
{
public:
    EdgeSpringForce(double springConstant, const std::vector<double>& edgeRestLengths)
        : m_springConstant(springConstant)
        , m_edgeRestLengths(edgeRestLengths) {
    }

    void apply(IFaceBody& body, double dt) override;

private:
    double m_springConstant;
    std::vector<double> m_edgeRestLengths;
};

/** Edge-aligned relative-velocity damping for each mesh edge */
class EdgeDampingForce : public IForceGenerator
{
public:
    explicit EdgeDampingForce(double dampingCoeff)
        : m_dampingCoeff(dampingCoeff) {
    }

    void apply(IFaceBody& body, double dt) override;

private:
    double m_dampingCoeff;
};


