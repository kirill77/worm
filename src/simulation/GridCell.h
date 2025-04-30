#pragma once

#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include "Molecule.h"
#include "MRNA.h"

// A single cell in the 3D grid representing the simulation space
class GridCell 
{
public:
    // Minimum possible resource level (to check with assertions)
    static constexpr double MIN_RESOURCE_LEVEL = 0.0;

    std::unordered_map<std::string, MPopulation> m_proteins;
    std::vector<std::shared_ptr<MRNA>> m_pMRNAs;
    
    // Constructor
    GridCell();
    
    // Helper to get or create protein population
    MPopulation& getOrCreateProtein(const std::string& sProteinName);
}; 