#pragma once

#include <vector>
#include <memory>
#include "BodyInterfaces.h"

// Semi-implicit Euler integrator operating on generic bodies
class PhysicsIntegrator
{
public:
    PhysicsIntegrator() = default;

    // Add a body to be integrated
    void addBody(std::shared_ptr<IFaceBody> body);

    // Integrate all registered bodies
    void step(double dt);

private:
    std::vector<std::shared_ptr<IFaceBody>> m_bodies;
};


