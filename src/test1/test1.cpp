#include "simulation/World.h"
#include "worm/Worm.h"
#include "log/ILog.h"
#include "simulation/ProteinWiki.h"
#include "visualization/Window.h"
#include "visualization/GPUWorld.h"
#include <memory>

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
    
    bool allTestsPassed = true;
    constexpr float fDtSec = 0.1f;  // 0.1 seconds per timestep
    float fCurrentTimeSec = 0.0f;   // Current simulation time in seconds
    
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
        pGPUWorld->drawMeshesIntoWindow();
    }

    if (allTestsPassed) {
        LOG_INFO("All development validation checks passed!");
    }

    return allTestsPassed ? 0 : 1;
}
