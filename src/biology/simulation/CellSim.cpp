#include "CellSim.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"
#include "geometry/geomHelpers/BVHMesh.h"
#include "chemistry/StringDict.h"
#include "physics/tensionSphere/TensionSphere.h"

CellSim::CellSim(std::shared_ptr<Cell> pCell)
    : m_pCell(pCell)
{
}

void CellSim::update(double dt)
{
    // Create tension sphere if it doesn't exist
    if (!m_pTensionSphere) {
        double volume = m_pCell->getInternalMedium().getVolumeMicroM();
        m_pTensionSphere = std::make_shared<TensionSphere>(2, volume);
    }
    
    // Update tension sphere physics
    m_pTensionSphere->makeTimeStep(dt);
    
    // Create cortex BVH if it doesn't exist
    if (!m_pCortexBVH)
    {
        auto pCortexMesh = m_pTensionSphere->getEdgeMesh();
        if (pCortexMesh)
        {
            m_pCortexBVH = std::make_shared<BVHMesh>(pCortexMesh);
            
            // Set the BVH mesh in the cortex organelle
            auto pCortexOrganelle = m_pCell->getOrganelle(StringDict::ID::ORGANELLE_CORTEX);
            if (pCortexOrganelle)
            {
                auto pCortex = std::dynamic_pointer_cast<Cortex>(pCortexOrganelle);
                if (pCortex)
                {
                    pCortex->setBVHMesh(m_pCortexBVH);
                }
            }
        }
    }

    m_pCell->update(dt);
}
