#pragma once

#include "biology/simulation/Organism.h"
#include "chemistry/DNA.h"
#include "biology/organelles/Chromosome.h"
#include "utils/dataCollector/DataCollector.h"
#include <memory>

class Worm : public Organism
{
private:
    std::shared_ptr<class Medium> createZygoteMedium();
    std::vector<Chromosome> initializeGenes();

    void setupDataCollector();  // Method to set up the data collector
    std::unique_ptr<DataCollector> m_pDataCollector;  // Data collector for simulation data
    float m_fTotalTime = 0.0f;  // Total simulation time in seconds

    // Validation thresholds based on experimental data
    static constexpr double ANTERIOR_POSTERIOR_RATIO_THRESHOLD = 3.0;  // Minimum ratio for proper PAR polarization
    static constexpr double NUCLEAR_SIZE_THRESHOLD = 0.8;             // Relative to initial size
    static constexpr double ASYMMETRIC_DIVISION_RATIO = 0.6;         // Ratio of anterior to posterior cell size

    // Development timing constants (in seconds)
    static constexpr float POLARITY_ESTABLISHMENT_END_SEC = 360.0f;    // 6 minutes
    static constexpr float POLARITY_MAINTENANCE_END_SEC = 600.0f;      // 10 minutes
    static constexpr float NUCLEAR_ENVELOPE_BREAKDOWN_SEC = 750.0f;    // 12.5 minutes
    static constexpr float SPINDLE_ASSEMBLY_START_SEC = 900.0f;        // 15 minutes
    static constexpr float DIVISION_START_SEC = 1100.0f;               // 18.3 minutes

public:
    Worm();

    // Override the simulateStep method from the base class
    void simulateStep(double dt) override;

    // Validation functions (time in seconds)
    bool validatePARPolarization(float fTimeSec) const;
    bool validateCellCycle(float fTimeSec) const;
    bool validateAsymmetricDivision(float fTimeSec) const;
    bool validateCentrosomeBehavior(float fTimeSec) const;
};

