#include "VisObjectContext.h"
#include "ConnectedMeshVis.h"
#include "biology/Organelle.h"
#include "biology/Organism.h"
#include "biology/Cell.h"
#include "biology/Cortex.h"
#include "visualization/gpu/GPUQueue.h"
#include "chemistry/StringDict.h"

static std::shared_ptr<ConnectedMeshVis> createCortexVis(
    std::shared_ptr<Organelle> pOrganelle,
    GPUQueue* pQueue)
{
    // get a connected mesh that shows how the cortex looks like
    auto pCortex = std::dynamic_pointer_cast<Cortex>(pOrganelle);
    auto pConnectedMesh = pCortex->getTensionSphere().getConnectedMesh();

    std::shared_ptr<ConnectedMeshVis> pCortexVis = std::make_shared<ConnectedMeshVis>(pQueue);
    pCortexVis->setConnectedMesh(pConnectedMesh);
    return pCortexVis;
}

std::shared_ptr<VisObjectContext> VisObjectContext::createForOrganelle(
    std::shared_ptr<Organelle> pOrganelle,
    StringDict::ID organelleId,
    GPUQueue* pQueue)
{
    auto pVisContext = std::make_shared<VisObjectContext>();
    
    // For now, we handle cortex specifically, but this can be extended
    // to handle other organelle types
    if (organelleId == StringDict::ID::ORGANELLE_CORTEX)
    {
        pVisContext->m_pObject = createCortexVis(pOrganelle, pQueue);
    }
    // TODO: Add other organelle visualization creation logic here
    
    pOrganelle->setVisObjectContext(pVisContext);
    return pVisContext;
} 