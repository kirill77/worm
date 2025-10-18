#pragma once

#include "BodyInterfaces.h"

// Semi-implicit Euler integrator operating on generic bodies
struct PhysicsIntegrator
{
    static void step(IFaceBody& body, double dt);
};


