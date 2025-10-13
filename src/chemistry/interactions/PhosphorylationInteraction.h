#pragma once

#include "MoleculeInteraction.h"
#include "chemistry/molecules/GridCell.h"
#include "chemistry/molecules/StringDict.h"

/**
 * Represents a phosphorylation interaction where one protein
 * adds a phosphate group to another, typically inactivating it.
 */
class PhosphorylationInteraction : public MoleculeInteraction
{
public:
    struct Parameters {
        double removalRate;         // Rate at which the kinase phosphorylates target
        double saturationConstant;  // Saturation constant for Hill-type kinetics
    };
    
    // Constructor
    PhosphorylationInteraction(StringDict::ID kinaseId, 
                              StringDict::ID targetId,
                              StringDict::ID phosphorylatedId,
                              const Parameters& params);
    
    // Apply phosphorylation to proteins in the cell
    bool apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const override;
    
private:
    StringDict::ID m_kinaseId;         // ID of the kinase protein
    StringDict::ID m_targetId;         // ID of the target protein
    StringDict::ID m_phosphorylatedId; // ID of the phosphorylated protein
    double m_removalRate;
    double m_saturationConstant;
}; 