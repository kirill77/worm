#include "simulation/World.h"
#include "simulation/Cell.h"
#include "worm/Worm.h"
#include "log/ILog.h"
#include "simulation/ProteinWiki.h"
#include "visualization/Window.h"
#include "visualization/GPUWorld.h"
#include "visualization/GPUStats.h"
#include "visHelpers/connectedMeshVis.h"
#include <memory>

std::shared_ptr<ConnectedMeshVis> createCortexVis(std::shared_ptr<Worm> pWorm, std::shared_ptr<Window> pWindow)
{
    // get a connected mesh that shows how the cortex looks like
    auto pCell = pWorm->getCells()[0];
    auto pCortex = pCell->getCortex();
    auto pConnectedMesh = pCortex->getTensionSphere().getConnectedMesh();

    std::shared_ptr<ConnectedMeshVis> pCortexVis = std::make_shared<ConnectedMeshVis>(pWindow);
    pCortexVis->setConnectedMesh(pConnectedMesh);
    return pCortexVis;
}

int main()
{
    // Initialize protein interaction data
    ProteinWiki::Initialize();

    // Create simulation objects
    std::shared_ptr<Worm> pWorm = std::make_shared<Worm>();
    World world(pWorm);

    // Create visualization window
    std::shared_ptr<Window> pWindow = std::make_shared<Window>();
    if (!pWindow->createWindowDevicAndSwapChain("Worm Simulation")) {
        LOG_ERROR("Failed to create window");
        return 1;
    }

    // Create GPU world for visualization
    std::shared_ptr<GPUWorld> pGPUWorld = std::make_shared<GPUWorld>(pWindow);

    std::shared_ptr<ConnectedMeshVis> pCortexVis = createCortexVis(pWorm, pWindow);
    pGPUWorld->addMesh(pCortexVis->getGPUMesh());

    bool allTestsPassed = true;
    constexpr float fDtSec = 0.1f;  // 0.1 seconds per timestep
    float fCurrentTimeSec = 0.0f;   // Current simulation time in seconds
    
    GPUStats gpuStats(pWindow->getDevice(), pWindow->createOrGetGPUQueue());

    // Main simulation loop
    while (true)
    {
        // Process window messages
        pWindow->processMessages();
        if (pWindow->shouldExit())
            break;
        
        // Simulate one step
        world.simulateStep(fDtSec);
        fCurrentTimeSec += fDtSec;

        // Run validation checks every 10 seconds
        if ((static_cast<uint32_t>(fCurrentTimeSec / fDtSec) + 1) % 100 == 0) {
            bool parValid = pWorm->validatePARPolarization(fCurrentTimeSec);
            bool cycleValid = pWorm->validateCellCycle(fCurrentTimeSec);
            bool divisionValid = pWorm->validateAsymmetricDivision(fCurrentTimeSec);
            
            if (!(parValid && cycleValid && divisionValid)) {
                allTestsPassed = false;
                LOG_ERROR("Validation failed at %.2lf sec", fCurrentTimeSec);
                break;
            }
        }

        // Render the visualization
        pCortexVis->updateGPUMesh(*pWindow->createOrGetGPUQueue());

        gpuStats.begin();
        pGPUWorld->drawMeshesIntoWindow();
        gpuStats.end();
        std::string s = gpuStats.getStats();
        s = s;
    }

    if (allTestsPassed) {
        LOG_INFO("All development validation checks passed!");
    }

    return allTestsPassed ? 0 : 1;
}
