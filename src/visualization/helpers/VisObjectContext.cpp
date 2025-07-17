#include "VisObjectContext.h"
#include "CortexVis.h"
#include "CentrosomeVis.h"
#include "biology/organelles/Organelle.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"
#include "biology/organelles/Centrosome.h"
#include "visualization/gpu/GPUQueue.h"
#include "chemistry/StringDict.h"

static std::shared_ptr<CortexVis> createCortexVis(
    std::shared_ptr<Organelle> pOrganelle,
    GPUQueue* pQueue)
{
    // get a connected mesh that shows how the cortex looks like
    auto pCortex = std::dynamic_pointer_cast<Cortex>(pOrganelle);
    std::shared_ptr<CortexVis> pCortexVis = std::make_shared<CortexVis>(pCortex, pQueue);
    return pCortexVis;
}

static std::shared_ptr<CentrosomeVis> createCentrosomeVis(
    std::shared_ptr<Organelle> pOrganelle,
    GPUQueue* pQueue)
{
    // create visualization for the centrosome
    auto pCentrosome = std::dynamic_pointer_cast<Centrosome>(pOrganelle);
    std::shared_ptr<CentrosomeVis> pCentrosomeVis = std::make_shared<CentrosomeVis>(pCentrosome, pQueue);
    return pCentrosomeVis;
}

std::shared_ptr<VisObjectContext> VisObjectContext::createForOrganelle(
    std::shared_ptr<Organelle> pOrganelle,
    StringDict::ID organelleId,
    GPUQueue* pQueue)
{
    auto pVisContext = std::make_shared<VisObjectContext>();
    
    // Handle different organelle types for visualization
    switch (organelleId)
    {
        case StringDict::ID::ORGANELLE_CORTEX:
            pVisContext->m_pObject = createCortexVis(pOrganelle, pQueue);
            break;
        case StringDict::ID::ORGANELLE_CENTROSOME:
            pVisContext->m_pObject = createCentrosomeVis(pOrganelle, pQueue);
            break;
        default:
            // TODO: Add other organelle visualization creation logic here
            break;
    }
    
    pOrganelle->setVisObjectContext(pVisContext);
    return pVisContext;
} 