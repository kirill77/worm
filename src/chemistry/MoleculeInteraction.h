#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include "Molecule.h"

// Forward declarations
class GridCell;
class ResourceDistributor;

/**
 * Base class for molecule interactions.
 * Simply performs actions on molecule populations without
 * describing the nature of the interaction.
 */
class MoleculeInteraction
{
public:
    // Molecular mechanisms (for informational purposes only)
    enum class Mechanism {
        PHOSPHORYLATION,    // Adding phosphate group
        DEPHOSPHORYLATION,  // Removing phosphate group
        BINDING,            // Physical binding
        CORTICAL_EXCLUSION, // Competitive binding to cortex
        RECRUITMENT,        // Recruiting to location
        DEGRADATION,        // Molecule degradation
        ACTIVATION,         // Conformational change activation
        INHIBITION          // Conformational change inhibition
    };
    
    // Basic constructor
    MoleculeInteraction(Mechanism mechanism, double atpCost) 
        : m_mechanism(mechanism)
        , m_atpCost(atpCost)
    {}
    
    // Virtual destructor for proper cleanup
    virtual ~MoleculeInteraction() = default;
    
    // Get mechanism (informational only)
    Mechanism getMechanism() const { return m_mechanism; }
    
    /**
     * Apply the interaction directly to the molecules in the cell
     * 
     * @param cell The grid cell containing molecules to act on
     * @param dt Time step in seconds
     * @param resDistributor Object to handle resource distribution
     * @return true if any changes were made, false otherwise
     */
    virtual bool apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const = 0;
    
protected:
    Mechanism m_mechanism;
    double m_atpCost;
};
