#pragma once

#include "Organism.h"
#include "DNA.h"
#include "Chromosome.h"

class Worm : public Organism
{
private:
    std::shared_ptr<class Medium> createZygoteMedium();
    std::vector<Chromosome> initializeGenes();

public:
    Worm();
};

