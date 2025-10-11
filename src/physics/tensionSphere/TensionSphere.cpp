#include "TensionSphere.h"
#include <cmath>
#include <algorithm>
#include <cassert>
#include "SoftBodyMeshAdapter.h"
#include "PhysicsIntegrator.h"
#include "VolumeConstraint.h"
#include "utils/log/ILog.h"

constexpr double PI = 3.14159265358979323846;

// TensionSphere implementation
TensionSphere::TensionSphere(uint32_t subdivisionLevel, double volume)
{
    // Initialize volume to the specified value
    m_fVolume = volume;
    
    // Create the base mesh with icosahedron and subdivisions using radius matching target volume
    double baseRadius = 1.0;
    if (m_fVolume > 0.0)
        baseRadius = std::cbrt(m_fVolume * 3.0 / (4.0 * PI));
    m_pMesh = std::make_shared<EdgeMesh>(baseRadius, subdivisionLevel);

    // Initialize physics simulation data
    initializePhysics();

    // Install default force generators (spring + damping) using current parameters
    m_forceGenerators.emplace_back(std::make_unique<EdgeSpringForce>(m_fSpringC));
    m_forceGenerators.emplace_back(std::make_unique<EdgeDampingForce>(m_fDampingCoeff));

    // Install volume constraint once (hard constraint by default: compliance=0)
    m_constraints.emplace_back(std::make_unique<VolumeConstraintXPBD>(m_fVolume, 0.0));
}

void TensionSphere::initializePhysics()
{
    const uint32_t vertexCount = m_pMesh->getVertexCount();
    
    // Initialize velocities to zero
    m_vertexVelocities.resize(vertexCount, double3(0, 0, 0));
    
    // Compute rest lengths for each edge from the mesh
    m_edgeRestLengths.clear();
    const uint32_t edgeCount = m_pMesh->getEdgeCount();
    m_edgeRestLengths.reserve(edgeCount);

    for (uint32_t e = 0; e < edgeCount; ++e)
    {
        auto edge = m_pMesh->getEdge(e);
        float3 pos1 = m_pMesh->getVertexPosition(edge.first);
        float3 pos2 = m_pMesh->getVertexPosition(edge.second);
        double restLength = length(pos2 - pos1);
        m_edgeRestLengths.push_back(restLength);
    }
}

void TensionSphere::makeTimeStep(double fDtSec)
{
    const uint32_t vertexCount = m_pMesh->getVertexCount();
    // Build adapter and force buffer
    SoftBodyMeshAdapter adapter(m_pMesh, m_vertexVelocities, m_edgeRestLengths);
    adapter.resizeForces();

    // Accumulate forces via generators
    for (auto& gen : m_forceGenerators)
        gen->apply(adapter, fDtSec);

    // Integrate motion using generic integrator
    PhysicsIntegrator::step(adapter, adapter.forces(), fDtSec);

    // Save pre-projection positions for velocity correction
    std::vector<double3> preProject(vertexCount);
    for (uint32_t i = 0; i < vertexCount; ++i)
        preProject[i] = double3(m_pMesh->getVertexPosition(i));

    // Apply registered constraints (projection)
    for (auto& c : m_constraints)
        c->project(adapter, fDtSec);

    // Post-projection velocity correction: v += (x_proj - x_pre) / dt
    if (fDtSec > 0.0)
    {
        for (uint32_t i = 0; i < vertexCount; ++i)
        {
            double3 xProj = double3(m_pMesh->getVertexPosition(i));
            m_vertexVelocities[i] += (xProj - preProject[i]) / fDtSec;
        }
    }
}

double TensionSphere::getVolume() const
{
    return m_fVolume;
}

void TensionSphere::setVolume(double volume)
{
    m_fVolume = volume;
}
