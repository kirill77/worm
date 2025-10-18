#include "PhysicsIntegrator.h"
#include "geometry/mesh/EdgeMesh.h"
#include <algorithm>

void PhysicsIntegrator::step(IFaceBody& body, double dt)
{
    if (dt <= 0.0) return;
    const uint32_t n = body.m_pMesh->getVertexCount();

    for (uint32_t i = 0; i < n; ++i)
    {
        auto& node = body.getVertex(i);
        const double m = std::max(1e-12, node.m_fMass);
        const double3 a = node.m_vForce / m;
        node.m_vVelocity += a * dt;
        const double3 xOld = double3(body.m_pMesh->getVertexPosition(i));
        const double3 xNew = xOld + node.m_vVelocity * dt;
        body.m_pMesh->setVertexPosition(i, float3(xNew));
        node.m_vForce = double3(0, 0, 0);
    }
}

