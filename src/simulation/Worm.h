#pragma once

#include "Organism.h"
#include "DNA.h"

class Worm : public Organism
{
private:
    std::shared_ptr<DNA> m_pDNA;
    std::shared_ptr<class Medium> createZygoteMedium();
    void initializeGenes();

public:
    Worm();
    
    // DNA access
    std::shared_ptr<DNA> getDNA() const { return m_pDNA; }
};

