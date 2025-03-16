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
        double recoveryRate;        // Rate at which phosphorylated proteins recover
        double saturationConstant;  // Saturation constant for Hill-type kinetics
    };
    
    // Constructor
    PhosphorylationInteraction(const std::string& kinase, 
                              const std::string& target, 
                              const Parameters& params);
    
    // Apply phosphorylation to proteins in the cell
    bool apply(GridCell& cell, double dt, double& atpConsumed) const override;
    
    // Handle recovery effects in neighboring cells
    bool applyNeighborEffects(GridCell& cell, 
                            std::vector<std::reference_wrapper<GridCell>>& neighborCells, 
                            double dt) const override;
    
private:
    // Kinase properties
    double m_removalRate;
    double m_recoveryRate;
    double m_saturationConstant;
}; 