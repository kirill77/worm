#pragma once

#include "ProteinInteraction.h"
#include "GridCell.h"
#include "StringDict.h"

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
    
    // Constructor - target is the base protein ID
    DephosphorylationInteraction(StringDict::ID targetId, 
                                StringDict::ID phosphorylatedId,
                                const Parameters& params);
    
    // Apply dephosphorylation to proteins in the cell
    bool apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const override;
    
private:
    StringDict::ID m_targetId;         // ID of the target protein
    StringDict::ID m_phosphorylatedId; // ID of the phosphorylated protein
    double m_recoveryRate;
}; 