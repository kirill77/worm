#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include "Protein.h"

// Forward declaration
class GridCell;

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
    ProteinInteraction(const std::string& proteinA, 
                     const std::string& proteinB,
                     Mechanism mechanism,
                     double atpCost) 
        : m_proteinA(proteinA)
        , m_proteinB(proteinB)
        , m_mechanism(mechanism)
        , m_atpCost(atpCost)
    {}
    
    // Virtual destructor for proper cleanup
    virtual ~ProteinInteraction() = default;
    
    // Get proteins involved in this interaction
    const std::string& getProteinA() const { return m_proteinA; }
    const std::string& getProteinB() const { return m_proteinB; }
    
    // Get mechanism (informational only)
    Mechanism getMechanism() const { return m_mechanism; }
    
    // Get ATP cost
    double getATPCost() const { return m_atpCost; }
    
    /**
     * Apply the interaction directly to the proteins in the cell
     * 
     * @param cell The grid cell containing proteins to act on
     * @param dt Time step in seconds
     * @param atpConsumed Reference to track ATP consumed by this interaction
     * @return true if any changes were made, false otherwise
     */
    virtual bool apply(GridCell& cell, double dt, double& atpConsumed) const = 0;
    
protected:
    std::string m_proteinA;
    std::string m_proteinB;
    Mechanism m_mechanism;
    double m_atpCost;
}; 