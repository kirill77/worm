#include "PhysicsIntegrator.h"
#include "geometry/mesh/EdgeMesh.h"
#include <algorithm>

void PhysicsIntegrator::addBody(std::shared_ptr<IFaceBody> body)
{
    m_bodies.push_back(body);
}

void PhysicsIntegrator::step(double dt)
{
    if (dt <= 0.0) return;

    for (auto& pBody : m_bodies)
    {
        const uint32_t n = pBody->m_pMesh->getVertexCount();

        for (uint32_t i = 0; i < n; ++i)
        {
            auto& node = pBody->getVertex(i);
            const double m = std::max(1e-12, node.m_fMass);
            const double3 a = node.m_vForce / m;
            node.m_vVelocity += a * dt;
            const double3 xOld = double3(pBody->m_pMesh->getVertexPosition(i));
            const double3 xNew = xOld + node.m_vVelocity * dt;
            pBody->m_pMesh->setVertexPosition(i, float3(xNew));
            node.m_vForce = double3(0, 0, 0);
        }
    }
}

