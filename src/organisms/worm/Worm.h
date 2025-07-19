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

