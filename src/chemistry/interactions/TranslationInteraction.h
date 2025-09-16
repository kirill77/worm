#pragma once

#include "MoleculeInteraction.h"
#include "chemistry/molecules/Molecule.h"
#include <vector>

/**
 * Handles translation of mRNA into proteins using tRNAs
 * 
 * This interaction represents the biological process where:
 * - mRNA provides the template
 * - tRNAs bring amino acids matching the codons
 * - Ribosomes facilitate the process
 * - Proteins are produced as a result
 */
class TranslationInteraction : public MoleculeInteraction
{
public:
    struct Parameters {
        double translationRate;     // Rate of protein production per mRNA per second
    };
    
    // Constructor taking mRNA molecule and translation parameters
    TranslationInteraction(const Molecule& mRNA, const Parameters& params);
    
    // Apply the translation interaction
    bool apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const override;
    
    // Get the mRNA being translated
    const Molecule& getMRNA() const { return m_mRNA; }

private:
    Molecule m_mRNA;                // The mRNA being translated
    double m_translationRate;       // Translation rate parameter
    
    // Helper methods
    void consumeTRNAs(GridCell& cell, 
                     const std::vector<std::pair<Molecule, uint32_t>>& requiredTRNAs,
                     double proteinAmount) const;
};
