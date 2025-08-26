#pragma once

#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include "Molecule.h"
#include "MRNA.h"
#include "TRNA.h"

// A single cell in the 3D grid representing the simulation space
class GridCell 
{
public:
    // Minimum possible resource level (to check with assertions)
    static constexpr double MIN_RESOURCE_LEVEL = 0.0;

    std::unordered_map<std::string, MPopulation> m_molecules;
    std::vector<std::shared_ptr<MRNA>> m_pMRNAs;
    std::vector<std::shared_ptr<TRNA>> m_pTRNAs;
    
    // Constructor
    GridCell();
    
    // Helper to get or create molecule population
    MPopulation& getOrCreateMolecule(const std::string& sMoleculeName);
    
    // mRNA management
    void updateMRNAs(double dt);  // Handle mRNA degradation and cleanup
    
    // tRNA management  
    void updateTRNAs(double dt);  // Handle tRNA charging and cleanup
}; 