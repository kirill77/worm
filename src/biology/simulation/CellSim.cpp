#include "CellSim.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"
#include "geometry/geomHelpers/BVHMesh.h"
#include "chemistry/StringDict.h"

CellSim::CellSim(std::shared_ptr<Cell> pCell)
    : m_pCell(pCell)
{
}

void CellSim::update(double dt)
{
    m_pCell->update(dt);
    
    // Create cortex BVH if it doesn't exist
    if (!m_pCortexBVH)
    {
        auto pCortexOrganelle = m_pCell->getOrganelle(StringDict::ID::ORGANELLE_CORTEX);
        if (pCortexOrganelle)
        {
            auto pCortex = std::dynamic_pointer_cast<Cortex>(pCortexOrganelle);
            if (pCortex)
            {
                auto pCortexMesh = pCortex->getMesh();
                if (pCortexMesh)
                {
                    m_pCortexBVH = std::make_shared<BVHMesh>(pCortexMesh);
                }
            }
        }
    }
}
