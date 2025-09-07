#include "Molecule.h"
#include "StringDict.h"
#include "GeneWiki.h"
#include "TRNA.h"
#include <algorithm>
#include <cmath>


std::shared_ptr<MPopulation> Molecule::translate(double dt, double moleculeAmount, double translationRate, 
                                               const std::unordered_map<Molecule, Population>& availableMolecules) const
{
    // Assert that this method is only called on RNA molecules
    assert(m_type == ChemicalType::MRNA && "translate() can only be called on mRNA molecules");
    
    // Check if we have enough RNA to produce protein
    if (moleculeAmount < 0.1) return nullptr;  // Threshold for translation

    // Calculate protein production based on translation rate and available RNA
    double fProteinAmount = translationRate * dt * moleculeAmount;

    // Get sequence from GeneWiki
    const std::string& sequence = GeneWiki::getInstance().getSequence(getName());

    // Check for required charged tRNAs (simplified - in reality would check entire sequence)
    bool hasRequiredTRNAs = true;
    std::vector<std::pair<const Molecule*, Population*>> usedTRNAs; // Keep track of used tRNAs for discharge
    
    for (size_t i = 0; i < sequence.length(); i += 3)
    {
        if (i + 2 >= sequence.length()) break;  // Incomplete codon
        
        std::string codon = sequence.substr(i, 3);
        bool codonMatched = false;

        // Convert codon to anticodon and find matching charged tRNAs
        std::string anticodon = TRNA::codonToAnticodon(codon);
        if (!anticodon.empty()) {
            // Get charged tRNA IDs that have this anticodon
            auto matchingTRNAIds = TRNA::getChargedTRNAsWithAnticodon(anticodon);
            for (StringDict::ID tRNAId : matchingTRNAIds) {
                Molecule tRNAMolecule(tRNAId, ChemicalType::TRNA);
                auto it = availableMolecules.find(tRNAMolecule);
                if (it != availableMolecules.end() && it->second.m_fNumber > 0.1) {
                    codonMatched = true;
                    usedTRNAs.emplace_back(&it->first, const_cast<Population*>(&it->second));
                    break;
                }
            }
        }

        if (!codonMatched)
        {
            hasRequiredTRNAs = false;
            break;
        }
    }

    if (!hasRequiredTRNAs)
    {
        return nullptr;
    }

    // Create new protein
    auto pProtein = std::make_shared<MPopulation>(Molecule(m_id, ChemicalType::PROTEIN), fProteinAmount);

    // Discharge used tRNAs by reducing their population and creating uncharged variants
    // Note: This is simplified - in reality we'd need to manage the charge state transitions
    // For now, we just reduce the charged tRNA population
    for (const auto& usedTRNA : usedTRNAs)
    {
        if (usedTRNA.second->m_fNumber > 0.1)
        {
            usedTRNA.second->m_fNumber -= 0.1; // Consume small amount for translation
        }
    }

    return pProtein;
}