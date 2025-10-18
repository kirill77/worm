#include "PhysicsIntegrator.h"
#include <algorithm>

void PhysicsIntegrator::step(IBody& body, std::vector<double3>& forces, double dt)
{
    if (dt <= 0.0) return;
    auto& N = body.nodes();
    const size_t n = N.size();
    if (forces.size() < n) forces.resize(n, double3(0,0,0));

    for (uint32_t i = 0; i < n; ++i)
    {
        const double m = std::max(1e-12, N.getMass(i));
        const double3 a = forces[i] / m;
        const double3 v = N.getVelocity(i) + a * dt;
        N.setVelocity(i, v);
        const double3 x = N.getPosition(i) + v * dt;
        N.setPosition(i, x);
        forces[i] = double3(0,0,0);
    }
}

