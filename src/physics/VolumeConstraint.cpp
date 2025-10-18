#include "VolumeConstraint.h"
#include "geometry/mesh/EdgeMesh.h"

double VolumeConstraintXPBD::computeSignedVolume() const
{
    const uint32_t faceCount = m_body.m_pMesh->getTriangleCount();
    if (faceCount == 0) return 0.0;
    double V = 0.0;
    for (uint32_t f = 0; f < faceCount; ++f)
    {
        uint3 tri = m_body.m_pMesh->getTriangleVertices(f);
        double3 a = double3(m_body.m_pMesh->getVertexPosition(tri.x));
        double3 b = double3(m_body.m_pMesh->getVertexPosition(tri.y));
        double3 c = double3(m_body.m_pMesh->getVertexPosition(tri.z));
        V += (1.0/6.0) * dot(a, cross(b, c));
    }
    return V;
}

void VolumeConstraintXPBD::project(double dt)
{
    const uint32_t faceCount = m_body.m_pMesh->getTriangleCount();
    if (faceCount == 0 || dt <= 0.0)
        return;

    // Accumulate per-vertex gradients dV/dx
    const uint32_t n = m_body.m_pMesh->getVertexCount();
    std::vector<double3> grad(n, double3(0,0,0));

    for (uint32_t f = 0; f < faceCount; ++f)
    {
        uint3 tri = m_body.m_pMesh->getTriangleVertices(f);
        uint32_t ia = tri.x, ib = tri.y, ic = tri.z;
        double3 a = double3(m_body.m_pMesh->getVertexPosition(ia));
        double3 b = double3(m_body.m_pMesh->getVertexPosition(ib));
        double3 c = double3(m_body.m_pMesh->getVertexPosition(ic));
        // dV/da = (1/6) (b x c), etc.
        grad[ia] += (1.0/6.0) * cross(b, c);
        grad[ib] += (1.0/6.0) * cross(c, a);
        grad[ic] += (1.0/6.0) * cross(a, b);
    }

    // Compute constraint value C(x) = V(x) - V_target
    const double V = computeSignedVolume();
    const double C = V - m_targetVolume;

    // Compute denominator sum w_i |grad_i|^2
    double denom = 0.0;
    for (uint32_t i = 0; i < n; ++i)
    {
        const double wi = 1.0 / std::max(1e-12, m_body.getVertex(i).m_fMass);
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
        const double wi = 1.0 / std::max(1e-12, m_body.getVertex(i).m_fMass);
        const double3 dx = -wi * deltaLambda * grad[i];
        const double3 xOld = double3(m_body.m_pMesh->getVertexPosition(i));
        const double3 xNew = xOld + dx;
        m_body.m_pMesh->setVertexPosition(i, float3(xNew));
    }
}


