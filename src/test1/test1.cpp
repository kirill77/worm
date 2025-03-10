#include "simulation/World.h"
#include "worm/Worm.h"
#include "log/ILog.h"

int main()
{
    std::shared_ptr<Worm> pWorm = std::make_shared<Worm>();
    World world(pWorm);
    
    bool allTestsPassed = true;
    constexpr float fDtSec = 0.1f;  // 0.1 seconds per timestep
    float fCurrentTimeSec = 0.0f;   // Current simulation time in seconds
    
    for (uint32_t u = 0; u < 12000; ++u)  // 12000 steps = 20 minutes
    {
        world.simulateStep(fDtSec);
        fCurrentTimeSec += fDtSec;

        // Run validation checks every 10 seconds
        if ((u + 1) % 100 == 0) {
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
    
    if (allTestsPassed) {
        LOG_INFO("All development validation checks passed!");
    }
    
    return allTestsPassed ? 0 : 1;
}
