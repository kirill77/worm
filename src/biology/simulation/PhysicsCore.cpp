#include "PhysicsCore.h"
#include <cmath>
#include <algorithm>
#include <cassert>
#include "SoftBodyMeshAdapter.h"
#include "physics/PhysicsIntegrator.h"
#include "physics/VolumeConstraint.h"
#include "utils/log/ILog.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"

void PhysicsCore::initialize(std::shared_ptr<Cell> pCell)
{
    m_pCell = pCell;

    // Pull mesh from cell's cortex
    auto pCortex = std::dynamic_pointer_cast<Cortex>(m_pCell->getOrganelle(StringDict::ID::ORGANELLE_CORTEX));
    assert(pCortex && "Cell must have a cortex organelle");
    m_pCortexMesh = pCortex->getEdgeMesh();

    // Pull volume from cell's internal medium
    m_fVolume = m_pCell->getInternalMedium().getVolumeMicroM();

    initializePhysics();

    // Create adapter once (reused across timesteps to avoid repeated allocations)
    m_pMeshAdapter = std::make_shared<SoftBodyMeshAdapter>(m_pCortexMesh);

    m_forceGenerators.emplace_back(std::make_unique<EdgeSpringForce>(m_fSpringC, m_edgeRestLengths));
    m_forceGenerators.emplace_back(std::make_unique<EdgeDampingForce>(m_fDampingCoeff));

    m_constraints.emplace_back(std::make_unique<VolumeConstraintXPBD>(m_fVolume, 0.0));
}

void PhysicsCore::initializePhysics()
{
    m_edgeRestLengths.clear();
    const uint32_t edgeCount = m_pCortexMesh->getEdgeCount();
    m_edgeRestLengths.reserve(edgeCount);
    for (uint32_t e = 0; e < edgeCount; ++e)
    {
        auto edge = m_pCortexMesh->getEdge(e);
        float3 pos1 = m_pCortexMesh->getVertexPosition(edge.first);
        float3 pos2 = m_pCortexMesh->getVertexPosition(edge.second);
        double restLength = length(pos2 - pos1);
        m_edgeRestLengths.push_back(restLength);
    }
}

void PhysicsCore::makeTimeStep(double fDtSec)
{
    // Pull updated volume from cell's internal medium
    m_fVolume = m_pCell->getInternalMedium().getVolumeMicroM();

    const uint32_t vertexCount = m_pCortexMesh->getVertexCount();

    for (auto& gen : m_forceGenerators)
        gen->apply(*m_pMeshAdapter, fDtSec);

    PhysicsIntegrator::step(*m_pMeshAdapter, fDtSec);

    std::vector<double3> preProject(vertexCount);
    for (uint32_t i = 0; i < vertexCount; ++i)
        preProject[i] = double3(m_pCortexMesh->getVertexPosition(i));

    for (auto& c : m_constraints)
        c->project(*m_pMeshAdapter, fDtSec);

    if (fDtSec > 0.0)
    {
        for (uint32_t i = 0; i < vertexCount; ++i)
        {
            double3 xProj = double3(m_pCortexMesh->getVertexPosition(i));
            m_pMeshAdapter->getVertex(i).m_vVelocity += (xProj - preProject[i]) / fDtSec;
        }
    }

    // Push updated mesh back to cortex
    auto pCortex = std::dynamic_pointer_cast<Cortex>(m_pCell->getOrganelle(StringDict::ID::ORGANELLE_CORTEX));
    pCortex->setMesh(m_pCortexMesh);
}


