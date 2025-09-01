#pragma once

#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include "Molecule.h"

#include "TRNA.h"

// A single cell in the 3D grid representing the simulation space
class GridCell 
{
public:
    // Minimum possible resource level (to check with assertions)
    static constexpr double MIN_RESOURCE_LEVEL = 0.0;

    std::unordered_map<Molecule, Population> m_molecules;
    std::vector<std::shared_ptr<TRNA>> m_pTRNAs;
    
    // Constructor
    GridCell();
    
    // Helper to get or create molecule population
    Population& getOrCreateMolPop(const Molecule& molecule);
    
    // Check if RNA molecules exist
    bool hasRNAs() const;
    
    // RNA management
    void updateRNAs(double dt);  // Handle RNA degradation and cleanup
    
    // tRNA management  
    void updateTRNAs(double dt);  // Handle tRNA charging and cleanup

private:
    // No longer need separate RNA storage - use m_molecules with RNA type
}; 