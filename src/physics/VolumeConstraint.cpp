#include "VolumeConstraint.h"

double VolumeConstraintXPBD::computeSignedVolume(const IBody& body)
{
    const auto& F = body.faces();
    const auto& N = body.nodes();
    const size_t faceCount = F.faceCount();
    if (faceCount == 0) return 0.0;
    double V = 0.0;
    for (uint32_t f = 0; f < faceCount; ++f)
    {
        uint3 tri = F.face(f);
        double3 a = N.getPosition(tri.x);
        double3 b = N.getPosition(tri.y);
        double3 c = N.getPosition(tri.z);
        V += (1.0/6.0) * dot(a, cross(b, c));
    }
    return V;
}

void VolumeConstraintXPBD::project(IBody& body, double dt)
{
    const auto& F = body.faces();
    auto& N = body.nodes();
    const size_t faceCount = F.faceCount();
    if (faceCount == 0 || dt <= 0.0)
        return;

    // Accumulate per-vertex gradients dV/dx
    const size_t n = N.size();
    std::vector<double3> grad(n, double3(0,0,0));

    for (uint32_t f = 0; f < faceCount; ++f)
    {
        uint3 tri = F.face(f);
        uint32_t ia = tri.x, ib = tri.y, ic = tri.z;
        double3 a = N.getPosition(ia);
        double3 b = N.getPosition(ib);
        double3 c = N.getPosition(ic);
        // dV/da = (1/6) (b x c), etc.
        grad[ia] += (1.0/6.0) * cross(b, c);
        grad[ib] += (1.0/6.0) * cross(c, a);
        grad[ic] += (1.0/6.0) * cross(a, b);
    }

    // Compute constraint value C(x) = V(x) - V_target
    const double V = computeSignedVolume(body);
    const double C = V - m_targetVolume;

    // Compute denominator sum w_i |grad_i|^2
    double denom = 0.0;
    for (uint32_t i = 0; i < n; ++i)
    {
        const double wi = 1.0 / std::max(1e-12, N.getMass(i));
        denom += wi * dot(grad[i], grad[i]);
    }
    if (denom <= 1e-20)
        return;

    // XPBD update
    const double alphaTilde = m_compliance / (dt * dt);
    // Flip sign so the update reduces C = V - V_target
    const double deltaLambda = (C - alphaTilde * m_lambda) / (denom + alphaTilde);
    m_lambda += deltaLambda;

    // Apply position corrections and velocity update
    for (uint32_t i = 0; i < n; ++i)
    {
        const double wi = 1.0 / std::max(1e-12, N.getMass(i));
        const double3 dx = -wi * deltaLambda * grad[i];
        const double3 xNew = N.getPosition(i) + dx;
        N.setPosition(i, xNew);
    }
}


