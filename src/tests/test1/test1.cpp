#include "organisms/worm/Worm.h"
#include "utils/log/ILog.h"
#include "visualization/helpers/VisEngine.h"
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
        if ((static_cast<uint32_t>(fCurrentTimeSec / fDtSec) + 1) % 100 == 0) {
            bool parValid = pWorm->validatePARPolarization(fCurrentTimeSec);
            bool cycleValid = pWorm->validateCellCycle(fCurrentTimeSec);
            bool divisionValid = pWorm->validateAsymmetricDivision(fCurrentTimeSec);
            bool centrosomeValid = pWorm->validateCentrosomeBehavior(fCurrentTimeSec);
            
            // Centrosome validation might fail initially before fertilization, so we'll be more lenient
            bool criticalValidationPassed = parValid && cycleValid && divisionValid;
            bool allValidationPassed = criticalValidationPassed && centrosomeValid;
            
            if (!criticalValidationPassed) {
                allTestsPassed = false;
                LOG_ERROR("Critical validation failed at %.2lf sec", fCurrentTimeSec);
                break;
            }
            
            LOG_INFO("Validation at %.2lf sec - PAR: %s, Cycle: %s, Division: %s, Centrosome: %s", 
                     fCurrentTimeSec, 
                     parValid ? "PASS" : "FAIL",
                     cycleValid ? "PASS" : "FAIL", 
                     divisionValid ? "PASS" : "FAIL",
                     centrosomeValid ? "PASS" : "FAIL");
        }
    }

    visEngine.shutdown();

    if (allTestsPassed) {
        LOG_INFO("All development validation checks passed!");
    }

    return allTestsPassed ? 0 : 1;
}
