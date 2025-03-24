#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include "Protein.h"

// Forward declarations
class GridCell;
class ResourceDistributor;

/**
 * Base class for protein interactions.
 * Simply performs actions on protein populations without
 * describing the nature of the interaction.
 */
class ProteinInteraction
{
public:
    // Molecular mechanisms (for informational purposes only)
    enum class Mechanism {
        PHOSPHORYLATION,    // Adding phosphate group
        DEPHOSPHORYLATION,  // Removing phosphate group
        BINDING,            // Physical binding
        CORTICAL_EXCLUSION, // Competitive binding to cortex
        RECRUITMENT,        // Recruiting to location
        DEGRADATION,        // Protein degradation
        ACTIVATION,         // Conformational change activation
        INHIBITION          // Conformational change inhibition
    };
    
    // Basic constructor
    ProteinInteraction(Mechanism mechanism, double atpCost) 
        : m_mechanism(mechanism)
        , m_atpCost(atpCost)
    {}
    
    // Virtual destructor for proper cleanup
    virtual ~ProteinInteraction() = default;
    
    // Get mechanism (informational only)
    Mechanism getMechanism() const { return m_mechanism; }
    
    // Get ATP cost
    double getATPCost() const { return m_atpCost; }
    
    /**
     * Apply the interaction directly to the proteins in the cell
     * 
     * @param cell The grid cell containing proteins to act on
     * @param dt Time step in seconds
     * @param resDistributor Object to handle resource distribution
     * @return true if any changes were made, false otherwise
     */
    virtual bool apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const = 0;
    
protected:
    Mechanism m_mechanism;
    double m_atpCost;
}; 