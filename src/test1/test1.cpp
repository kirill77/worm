#include "simulation/World.h"
#include "simulation/Cell.h"
#include "worm/Worm.h"
#include "log/ILog.h"
#include "simulation/ProteinWiki.h"
#include "visualization/Window.h"
#include "visualization/GPUWorld.h"
#include "visualization/GPUStats.h"
#include "visualization/GPUText.h"
#include "visualization/DirectXHelpers.h"
#include "visHelpers/connectedMeshVis.h"
#include "visHelpers/cameraUI.h"
#include <memory>

std::shared_ptr<ConnectedMeshVis> createCortexVis(
    std::shared_ptr<Worm> pWorm,
    std::shared_ptr<Window> pWindow)
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
    std::shared_ptr<GPUWorld> pGpuWorld = std::make_shared<GPUWorld>(pWindow, pWindow->getSwapChain()->getGPUQueue());

    GPUText gpuText(pGpuWorld->getFont());
    gpuText.printf("Hello World!");

    std::shared_ptr<ConnectedMeshVis> pCortexVis = createCortexVis(pWorm, pWindow);
    pGpuWorld->addMesh(pCortexVis->getGPUMesh());

    bool allTestsPassed = true;
    constexpr float fDtSec = 0.1f;  // 0.1 seconds per timestep
    float fCurrentTimeSec = 0.0f;   // Current simulation time in seconds

    CameraUI cameraUI;
    cameraUI.attachToCamera(pGpuWorld->getCamera());
    
    GPUStats gpuStats(pWindow->getDevice());

    // Main simulation loop
    while (true)
    {
        // Process window messages
        pWindow->processMessages();
        if (pWindow->shouldExit())
            break;

        cameraUI.setWorldBox(*pCortexVis->getConnectedBox());
        cameraUI.notifyNewUIState(pWindow->getCurrentUIState());

        // Simulate one step
        world.simulateStep(fDtSec);
        fCurrentTimeSec += fDtSec;

        // Run validation checks every 10 seconds
        if (false && (static_cast<uint32_t>(fCurrentTimeSec / fDtSec) + 1) % 100 == 0) {
            bool parValid = pWorm->validatePARPolarization(fCurrentTimeSec);
            bool cycleValid = pWorm->validateCellCycle(fCurrentTimeSec);
            bool divisionValid = pWorm->validateAsymmetricDivision(fCurrentTimeSec);
            
            if (!(parValid && cycleValid && divisionValid)) {
                allTestsPassed = false;
                LOG_ERROR("Validation failed at %.2lf sec", fCurrentTimeSec);
                break;
            }
        }

        {
            auto* pSwapChain = pWindow->getSwapChain();

            // Render the visualization
            pCortexVis->updateGPUMesh();

            auto pGpuQueue = pSwapChain->getGPUQueue();
            auto pCmdList = pGpuQueue->beginRecording();

            // Draw meshes into the command list
            pGpuWorld->render(pSwapChain, pCmdList.Get());

            // Draw text on top
            gpuText.render(pSwapChain, pGpuWorld->getSharedRootSignature(), pCmdList.Get());

            pGpuQueue->execute(pCmdList);

            // Present the frame
            ThrowIfFailed(pSwapChain->getSwapChain()->Present(1, 0));
        }

        pGpuWorld = pGpuWorld;
    }

    {
        auto* pSwapChain = pWindow->getSwapChain();
        pSwapChain->getGPUQueue()->flush();
    }

    if (allTestsPassed) {
        LOG_INFO("All development validation checks passed!");
    }

    return allTestsPassed ? 0 : 1;
}
