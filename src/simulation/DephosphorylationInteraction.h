#pragma once

#include "ProteinInteraction.h"
#include "GridCell.h"

/**
 * Represents a dephosphorylation interaction where phosphorylated proteins
 * lose their phosphate group and return to their original state.
 */
class DephosphorylationInteraction : public ProteinInteraction
{
public:
    struct Parameters {
        double recoveryRate;        // Rate at which phosphorylated proteins recover
    };
    
    // Constructor - target is the base protein name (e.g. "PAR-2")
    DephosphorylationInteraction(const std::string& targetName, 
                                const Parameters& params);
    
    // Apply dephosphorylation to proteins in the cell
    bool apply(GridCell& cell, double dt, double& atpConsumed) const override;
    
private:
    std::string m_targetName;    // Name of the target protein
    double m_recoveryRate;
}; 