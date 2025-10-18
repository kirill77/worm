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

    // Create adapter once (reused across timesteps to avoid repeated allocations)
    m_pMeshAdapter = std::make_shared<SoftBodyMeshAdapter>(m_pCortexMesh);

    // Register body with integrator
    m_integrator.addBody(m_pMeshAdapter);

    m_forceGenerators.emplace_back(std::make_unique<EdgeSpringForce>(*m_pMeshAdapter, m_fSpringC));
    m_forceGenerators.emplace_back(std::make_unique<EdgeDampingForce>(*m_pMeshAdapter, m_fDampingCoeff));

    m_constraints.emplace_back(std::make_unique<VolumeConstraintXPBD>(*m_pMeshAdapter, m_fVolume, 0.0));
}

void PhysicsCore::makeTimeStep(double fDtSec)
{
    // Pull updated volume from cell's internal medium
    m_fVolume = m_pCell->getInternalMedium().getVolumeMicroM();

    const uint32_t vertexCount = m_pCortexMesh->getVertexCount();

    for (auto& gen : m_forceGenerators)
        gen->apply(fDtSec);

    m_integrator.step(fDtSec);

    std::vector<double3> preProject(vertexCount);
    for (uint32_t i = 0; i < vertexCount; ++i)
        preProject[i] = double3(m_pCortexMesh->getVertexPosition(i));

    for (auto& c : m_constraints)
        c->project(fDtSec);

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


