#pragma once

#include "ProteinInteraction.h"
#include "GridCell.h"

/**
 * Represents a phosphorylation interaction where one protein
 * adds a phosphate group to another, typically inactivating it.
 */
class PhosphorylationInteraction : public ProteinInteraction
{
public:
    struct Parameters {
        double removalRate;         // Rate at which the kinase phosphorylates target
        double saturationConstant;  // Saturation constant for Hill-type kinetics
    };
    
    // Constructor
    PhosphorylationInteraction(const std::string& kinaseName, 
                              const std::string& targetName,
                              const Parameters& params);
    
    // Apply phosphorylation to proteins in the cell
    bool apply(GridCell& cell, double dt, ResourceAllocation& resDistributor) const override;
    
private:
    std::string m_kinaseName;    // Name of the kinase protein
    std::string m_targetName;    // Name of the target protein
    std::string m_phosphorylatedName; // Cached name of phosphorylated protein
    double m_removalRate;
    double m_saturationConstant;
}; 