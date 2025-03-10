#pragma once

#include "simulation/Organism.h"
#include "simulation/DNA.h"
#include "simulation/Chromosome.h"

class Worm : public Organism
{
private:
    std::shared_ptr<class Medium> createZygoteMedium();
    std::vector<Chromosome> initializeGenes();

    // Validation thresholds based on experimental data
    static constexpr double ANTERIOR_POSTERIOR_RATIO_THRESHOLD = 3.0;  // Minimum ratio for proper PAR polarization
    static constexpr double NUCLEAR_SIZE_THRESHOLD = 0.8;             // Relative to initial size
    static constexpr double ASYMMETRIC_DIVISION_RATIO = 0.6;         // Ratio of anterior to posterior cell size

public:
    Worm();

    // Validation functions
    bool validatePARPolarization(uint32_t timestep) const;
    bool validateCellCycle(uint32_t timestep) const;
    bool validateAsymmetricDivision(uint32_t timestep) const;
};

