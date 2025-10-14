#include "CellSim.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"
#include "biology/simulation/TensionSphere.h"

CellSim::CellSim(std::shared_ptr<Cell> pCell)
    : m_pCell(pCell)
{
    // Wire cortex mesh and pass it into TensionSphere
    auto pCortex = std::dynamic_pointer_cast<Cortex>(m_pCell->getOrganelle(StringDict::ID::ORGANELLE_CORTEX));

    // Cortex constructs its mesh in its constructor; reuse it directly
    auto pCortexMesh = pCortex->getEdgeMesh();
    m_pTensionSphere = std::make_shared<TensionSphere>(pCortexMesh, m_pCell->getInternalMedium().getVolumeMicroM());
}

void CellSim::update(const TimeContext& time)
{
    // Advance physics on the tension sphere first
    assert(m_pTensionSphere && m_pCell);
    double fVolume = m_pCell->getInternalMedium().getVolumeMicroM();
    m_pTensionSphere->setVolume(fVolume);
    m_pTensionSphere->makeTimeStep(time.m_deltaTSec);

    // Refresh cortex mesh BVH if needed
    auto pCortex = std::dynamic_pointer_cast<Cortex>(m_pCell->getOrganelle(StringDict::ID::ORGANELLE_CORTEX));
    pCortex->setMesh(m_pTensionSphere->getEdgeMesh());

    m_pCell->update(time.m_deltaTSec);
}
