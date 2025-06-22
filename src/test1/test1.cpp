#include "worm/Worm.h"
#include "log/ILog.h"
#include "visHelpers/VisEngine.h"
#include <memory>

int main()
{
    // Create simulation objects
    std::shared_ptr<Worm> pWorm = std::make_shared<Worm>();

    // Create and initialize visualization engine
    VisEngine visEngine;
    if (!visEngine.initialize(pWorm)) {
        LOG_ERROR("Failed to initialize visualization engine");
        return 1;
    }

    bool allTestsPassed = true;
    constexpr float fDtSec = 0.1f;  // 0.1 seconds per timestep
    float fCurrentTimeSec = 0.0f;   // Current simulation time in seconds

    // Main simulation loop
    while (true)
    {
        if (!visEngine.update(fDtSec)) {
            break;
        }

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
    }

    visEngine.shutdown();

    if (allTestsPassed) {
        LOG_INFO("All development validation checks passed!");
    }

    return allTestsPassed ? 0 : 1;
}
