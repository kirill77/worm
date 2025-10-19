#include "PhysicsIntegrator.h"
#include "geometry/mesh/EdgeMesh.h"
#include <algorithm>

void PhysicsIntegrator::addBody(std::shared_ptr<PhysicsMesh> body)
{
    m_bodies.push_back(body);
}

void PhysicsIntegrator::addForceGenerator(std::unique_ptr<IForceGenerator> generator)
{
    m_forceGenerators.push_back(std::move(generator));
}

void PhysicsIntegrator::addConstraint(std::shared_ptr<IConstraint> constraint)
{
    m_constraints.push_back(constraint);
}

void PhysicsIntegrator::step(double dt)
{
    if (dt <= 0.0) return;

    // Step 1: Apply all force generators
    for (auto& gen : m_forceGenerators)
        gen->apply(dt);

    // Step 2: Semi-implicit Euler integration
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
        }
    }

    // Step 3: Save pre-constraint positions for velocity correction
    std::vector<std::vector<double3>> preProjectPositions;
    preProjectPositions.reserve(m_bodies.size());
    for (auto& pBody : m_bodies)
    {
        const uint32_t n = pBody->m_pMesh->getVertexCount();
        std::vector<double3> positions(n);
        for (uint32_t i = 0; i < n; ++i)
            positions[i] = double3(pBody->m_pMesh->getVertexPosition(i));
        preProjectPositions.push_back(std::move(positions));
    }

    // Step 4: Apply all constraints (XPBD position corrections)
    for (auto& c : m_constraints)
        c->project(dt);

    // Step 5: Update velocities based on constraint-induced position changes
    for (size_t bodyIdx = 0; bodyIdx < m_bodies.size(); ++bodyIdx)
    {
        auto& pBody = m_bodies[bodyIdx];
        const uint32_t n = pBody->m_pMesh->getVertexCount();
        const auto& prePos = preProjectPositions[bodyIdx];

        for (uint32_t i = 0; i < n; ++i)
        {
            double3 xProj = double3(pBody->m_pMesh->getVertexPosition(i));
            pBody->getVertex(i).m_vVelocity += (xProj - prePos[i]) / dt;
        }
    }

    // Step 6: Clear forces for next timestep
    for (auto& pBody : m_bodies)
    {
        const uint32_t n = pBody->m_pMesh->getVertexCount();
        for (uint32_t i = 0; i < n; ++i)
            pBody->getVertex(i).m_vForce = double3(0, 0, 0);
    }
}

