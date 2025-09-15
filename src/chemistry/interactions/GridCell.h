#pragma once

#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include "chemistry/molecules/Molecule.h"

// A single cell in the 3D grid representing the simulation space
class GridCell 
{
public:
    // Minimum possible resource level (to check with assertions)
    static constexpr double MIN_RESOURCE_LEVEL = 0.0;

    std::unordered_map<Molecule, Population> m_molecules;
    
    // Constructor
    GridCell();
    
    // Helper to get or create molecule population
    Population& getOrCreateMolPop(const Molecule& molecule);
    
    // Check if mRNA molecules exist
    bool hasMRNAs() const;
    
    // mRNA management
    void updateMRNAs(double dt);  // Handle mRNA degradation and cleanup
    
    // tRNA management  
    void updateTRNAs(double dt);  // Handle tRNA charging transitions (uncharged -> charged)

private:
    // No longer need separate RNA storage - use m_molecules with MRNA type
}; 