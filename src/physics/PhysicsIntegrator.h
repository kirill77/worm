#pragma once

#include <vector>
#include <memory>
#include "PhysicsMesh.h"
#include "ForceGenerator.h"
#include "PhysicsConstraints.h"

// Physics integrator managing the complete simulation pipeline:
// force application, integration, constraint projection, and velocity correction
class PhysicsIntegrator
{
public:
    PhysicsIntegrator() = default;

    // Add a body to be integrated
    void addBody(std::shared_ptr<PhysicsMesh> body);

    // Add a force generator to the simulation
    void addForceGenerator(std::unique_ptr<IForceGenerator> generator);

    // Add a constraint to the simulation
    void addConstraint(std::shared_ptr<IConstraint> constraint);

    // Execute complete physics pipeline: forces -> integration -> constraints
    void step(double dt);

private:
    std::vector<std::shared_ptr<PhysicsMesh>> m_bodies;
    std::vector<std::unique_ptr<IForceGenerator>> m_forceGenerators;
    std::vector<std::shared_ptr<IConstraint>> m_constraints;
};


