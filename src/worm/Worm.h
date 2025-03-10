#pragma once

#include "simulation/Organism.h"
#include "simulation/DNA.h"
#include "simulation/Chromosome.h"

class Worm : public Organism
{
private:
    std::shared_ptr<class Medium> createZygoteMedium();
    std::vector<Chromosome> initializeGenes();

public:
    Worm();
};

