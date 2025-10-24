#include "PhysicsIntegrator.h"
#include "geometry/mesh/TriangleMesh.h"
#include <algorithm>

namespace {
    // Semi-implicit Euler integration for a single physics vertex
    inline void integratePhysVertex(PhysVertex& vertex, double dt) {
        const double m = std::max(1e-12, vertex.m_fMass);
        const double3 a = vertex.m_vForce / m;
        vertex.m_vVelocity += a * dt;
    }
    
    // Update position using integrated velocity
    inline void updatePosition(float3& position, const double3& velocity, double dt) {
        const double3 xOld = double3(position);
        const double3 xNew = xOld + velocity * dt;
        position = float3(xNew);
    }
}

void PhysicsIntegrator::addBody(std::shared_ptr<PhysicsMesh> body)
{
    m_bodies.push_back(body);
}

void PhysicsIntegrator::addCentrosome(std::shared_ptr<PhysCentrosome> centrosome)
{
    m_centrosomes.push_back(centrosome);
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
        gen->apply();

    // Step 2: Semi-implicit Euler integration for centrosomes
    for (auto& pCentrosome : m_centrosomes)
    {
        auto& physVertex = pCentrosome->getPhysVertex();
        integratePhysVertex(physVertex, dt);
        updatePosition(pCentrosome->getToNormalizedCell().m_translation, physVertex.m_vVelocity, dt);
    }

    // Step 3: Semi-implicit Euler integration for mesh bodies
    for (auto& pBody : m_bodies)
    {
        const uint32_t n = pBody->m_pMesh->getVertices()->getVertexCount();

        for (uint32_t i = 0; i < n; ++i)
        {
            auto& vertex = pBody->getVertex(i);
            integratePhysVertex(vertex, dt);
            
            // Update position in mesh geometry
            float3 pos = pBody->m_pMesh->getVertices()->getVertexPosition(i);
            updatePosition(pos, vertex.m_vVelocity, dt);
            pBody->m_pMesh->getVertices()->setVertexPosition(i, pos);
        }
    }

    // Step 4: Save pre-constraint positions for velocity correction
    std::vector<std::vector<double3>> preProjectPositions;
    preProjectPositions.reserve(m_bodies.size());
    for (auto& pBody : m_bodies)
    {
        const uint32_t n = pBody->m_pMesh->getVertices()->getVertexCount();
        std::vector<double3> positions(n);
        for (uint32_t i = 0; i < n; ++i)
            positions[i] = double3(pBody->m_pMesh->getVertices()->getVertexPosition(i));
        preProjectPositions.push_back(std::move(positions));
    }

    // Step 5: Apply all constraints (XPBD position corrections)
    for (auto& c : m_constraints)
        c->project(dt);

    // Step 6: Update velocities based on constraint-induced position changes
    for (size_t bodyIdx = 0; bodyIdx < m_bodies.size(); ++bodyIdx)
    {
        auto& pBody = m_bodies[bodyIdx];
        const uint32_t n = pBody->m_pMesh->getVertices()->getVertexCount();
        const auto& prePos = preProjectPositions[bodyIdx];

        for (uint32_t i = 0; i < n; ++i)
        {
            double3 xProj = double3(pBody->m_pMesh->getVertices()->getVertexPosition(i));
            pBody->getVertex(i).m_vVelocity += (xProj - prePos[i]) / dt;
        }
    }

    // Step 7: Clear forces for next timestep
    for (auto& pCentrosome : m_centrosomes)
    {
        pCentrosome->getPhysVertex().m_vForce = double3(0, 0, 0);
    }
    
    for (auto& pBody : m_bodies)
    {
        const uint32_t n = pBody->m_pMesh->getVertices()->getVertexCount();
        for (uint32_t i = 0; i < n; ++i)
            pBody->getVertex(i).m_vForce = double3(0, 0, 0);
    }
}

