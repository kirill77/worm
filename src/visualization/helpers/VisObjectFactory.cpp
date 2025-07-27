#include "VisObjectFactory.h"
#include "visualization/gpu/IVisObject.h"
#include "CortexVis.h"
#include "CentrosomeVis.h"
#include "biology/organelles/Organelle.h"
#include "biology/organelles/Cortex.h"
#include "biology/organelles/Centrosome.h"
#include "visualization/gpu/GPUQueue.h"
#include "chemistry/StringDict.h"

static std::shared_ptr<CortexVis> createCortexVis(
    std::shared_ptr<Organelle> pOrganelle,
    GPUQueue* pQueue)
{
    // Get a connected mesh that shows how the cortex looks like
    auto pCortex = std::dynamic_pointer_cast<Cortex>(pOrganelle);
    std::shared_ptr<CortexVis> pCortexVis = std::make_shared<CortexVis>(pCortex, pQueue);
    return pCortexVis;
}

static std::shared_ptr<CentrosomeVis> createCentrosomeVis(
    std::shared_ptr<Organelle> pOrganelle,
    GPUQueue* pQueue)
{
    // Create visualization for the centrosome
    auto pCentrosome = std::dynamic_pointer_cast<Centrosome>(pOrganelle);
    std::shared_ptr<CentrosomeVis> pCentrosomeVis = std::make_shared<CentrosomeVis>(pCentrosome, pQueue);
    return pCentrosomeVis;
}

std::shared_ptr<IVisObject> VisObjectFactory::createForOrganelle(
    std::shared_ptr<Organelle> pOrganelle,
    StringDict::ID organelleId,
    GPUQueue* pQueue)
{
    std::shared_ptr<IVisObject> pVisObject = nullptr;
    
    // Handle different organelle types for visualization
    switch (organelleId)
    {
        case StringDict::ID::ORGANELLE_CORTEX:
            pVisObject = createCortexVis(pOrganelle, pQueue);
            break;
        case StringDict::ID::ORGANELLE_CENTROSOME:
            pVisObject = createCentrosomeVis(pOrganelle, pQueue);
            break;
        default:
            // TODO: Add other organelle visualization creation logic here
            break;
    }
    
    // Set the visualization object in the organelle
    if (pVisObject)
    {
        pOrganelle->setVisObject(pVisObject);
    }
    
    return pVisObject;
} 