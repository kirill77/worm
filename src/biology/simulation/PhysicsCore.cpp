#include "PhysicsCore.h"
#include <cmath>
#include <algorithm>
#include <cassert>
#include "physics/PhysicsMesh.h"
#include "physics/PhysicsIntegrator.h"
#include "physics/VolumeConstraint.h"
#include "physics/DyneinPullingForce.h"
#include "physics/PhysCentrosome.h"
#include "chemistry/molecules/simConstants.h"
#include "utils/log/ILog.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"
#include "biology/organelles/Centrosome.h"

void PhysicsCore::initialize(std::shared_ptr<Cell> pCell)
{
    m_pCell = pCell;

    // Pull mesh from cell's cortex and create physics mesh
    auto pCortex = std::dynamic_pointer_cast<Cortex>(m_pCell->getOrganelle(StringDict::ID::ORGANELLE_CORTEX));
    assert(pCortex && "Cell must have a cortex organelle");
    auto pCortexMesh = pCortex->getTriangleMesh();
    m_pCortexAdapter = std::make_shared<PhysicsMesh>(pCortexMesh);

    // Register body with integrator
    m_integrator.addBody(m_pCortexAdapter);

    // Get centrosomes for dynein force
    auto pCentrosome = std::dynamic_pointer_cast<Centrosome>(m_pCell->getOrganelle(StringDict::ID::ORGANELLE_CENTROSOME));
    auto pPhysCentrosome = pCentrosome->getPhysCentrosome();
    m_centrosomes.push_back(pPhysCentrosome);

    // Register force generators with integrator
    m_integrator.addForceGenerator(std::make_unique<EdgeSpringForce>(*m_pCortexAdapter, m_fSpringC));
    m_integrator.addForceGenerator(std::make_unique<EdgeDampingForce>(*m_pCortexAdapter, m_fDampingCoeff));
    m_integrator.addForceGenerator(
        std::make_unique<DyneinPullingForce>(*m_pCortexAdapter, m_centrosomes, MoleculeConstants::DYNEIN_PULLING_FORCE_PICONEWTONS));

    // Pull volume from cell's internal medium
    double fVolume = m_pCell->getInternalMedium().getVolumeMicroM();
    // Register constraints with integrator
    m_pVolumeConstraint = std::make_shared<VolumeConstraintXPBD>(*m_pCortexAdapter, fVolume, 0.0);
    m_integrator.addConstraint(m_pVolumeConstraint);
}

void PhysicsCore::makeTimeStep(double fDtSec)
{
    // Update target volume from cell's internal medium
    double fVolume = m_pCell->getInternalMedium().getVolumeMicroM();
    m_pVolumeConstraint->setTargetVolume(fVolume);

    // Execute complete physics pipeline
    m_integrator.step(fDtSec);

    // Push updated mesh back to cortex
    auto pCortex = std::dynamic_pointer_cast<Cortex>(m_pCell->getOrganelle(StringDict::ID::ORGANELLE_CORTEX));
    pCortex->setTriangleMesh(m_pCortexAdapter->m_pMesh);
}


