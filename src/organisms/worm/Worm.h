#pragma once

#include "biology/simulation/Organism.h"
#include "biology/simulation/TimeContext.h"
#include "chemistry/molecules/DNA.h"
#include "biology/organelles/Chromosome.h"
#include "biology/dataCollector/DataCollector.h"
#include <memory>

class Worm : public Organism
{
private:
    std::shared_ptr<class Medium> createZygoteMedium();
    std::vector<Chromosome> initializeGenes();
    void addMaternalTRNAs(Medium& medium, const float3& position);

    void setupDataCollector();  // Method to set up the data collector
    std::unique_ptr<DataCollector> m_pDataCollector;  // Data collector for simulation data
    float m_fTotalTime = 0.0f;  // Total simulation time in seconds

public:
    Worm();

    // Override the simulateStep method from the base class
    void simulateStep(const TimeContext& time) override;

    // Validation functions (time in seconds)
    bool validatePARPolarization(float fTimeSec) const;
    bool validateCellCycle(float fTimeSec) const;
    bool validateAsymmetricDivision(float fTimeSec) const;
    bool validateCentrosomeBehavior(float fTimeSec) const;
    bool validateGammaTubulinLevels(float fTimeSec) const;
};

