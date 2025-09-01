#include "Molecule.h"
#include "StringDict.h"
#include "TRNA.h"
#include "GeneWiki.h"
#include <algorithm>
#include <cmath>

// Constructor with name and type - automatically tries to use ID optimization
Molecule::Molecule(const std::string& name, ChemicalType type) : m_type(type)
{
    // Try to convert string to StringDict ID first for optimization
    m_id = StringDict::stringToId(name);

    assert(m_type != ChemicalType::OTHER);

    if (m_id != StringDict::ID::eUNKNOWN) {
        // Found matching ID - use optimized storage (empty string)
        m_sName.clear();  // Explicitly clear for safety
    } else {
        // Unknown molecule - store the string
        m_sName = name;
    }
}

std::shared_ptr<MPopulation> Molecule::translate(double dt, double moleculeAmount, double translationRate, 
                                               const std::vector<std::shared_ptr<TRNA>>& availableTRNAs) const
{
    // Assert that this method is only called on RNA molecules
    assert(m_type == ChemicalType::MRNA && "translate() can only be called on mRNA molecules");
    
    // Check if we have enough RNA to produce protein
    if (moleculeAmount < 0.1) return nullptr;  // Threshold for translation

    // Calculate protein production based on translation rate and available RNA
    double fProteinAmount = translationRate * dt * moleculeAmount;

    // Get sequence from GeneWiki
    const std::string& sequence = GeneWiki::getInstance().getSequence(getName());

    // Check for required tRNAs (simplified - in reality would check entire sequence)
    bool hasRequiredTRNAs = true;
    for (size_t i = 0; i < sequence.length(); i += 3)
    {
        if (i + 2 >= sequence.length()) break;  // Incomplete codon
        
        std::string codon = sequence.substr(i, 3);  // substr returns a new string
        bool codonMatched = false;

        for (const auto& tRNA : availableTRNAs)
        {
            if (tRNA->isCharged() && tRNA->matchesCodon(codon))
            {
                codonMatched = true;
                break;
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
    auto pProtein = std::make_shared<MPopulation>(Molecule(getName(), ChemicalType::PROTEIN), fProteinAmount);

    // Discharge used tRNAs (simplified)
    for (const auto& tRNA : availableTRNAs)
    {
        if (tRNA->isCharged())
        {
            const_cast<TRNA*>(tRNA.get())->discharge();
        }
    }

    return pProtein;
}