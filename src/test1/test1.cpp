#include "simulation/World.h"
#include "worm/Worm.h"
#include "log/ILog.h"

int main()
{
    std::shared_ptr<Worm> pWorm = std::make_shared<Worm>();
    World world(pWorm);
    
    bool allTestsPassed = true;
    
    for (uint32_t u = 0; u < 1000; ++u)
    {
        world.simulateStep(0.1);
        
        // Run validation checks every 50 timesteps
        if (u % 50 == 0) {
            bool parValid = pWorm->validatePARPolarization(u);
            bool cycleValid = pWorm->validateCellCycle(u);
            bool divisionValid = pWorm->validateAsymmetricDivision(u);
            
            if (!(parValid && cycleValid && divisionValid)) {
                allTestsPassed = false;
                LOG_ERROR("Validation failed at timestep");
                break;
            }
        }
    }
    
    if (allTestsPassed) {
        LOG_INFO("All development validation checks passed!");
    }
        
    return allTestsPassed ? 0 : 1;
}
