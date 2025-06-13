#include "pch.h"
#include "MRNA.h"
#include "TRNA.h"
#include "Molecule.h"
#include "GeneWiki.h"
#include <algorithm>

std::shared_ptr<MPopulation> MRNA::translate(double dt, const std::vector<std::shared_ptr<TRNA>>& availableTRNAs) const
{
    // Check if we have enough mRNA to produce protein
    if (m_fNumber < 0.1) return nullptr;  // Threshold for translation

    // Calculate protein production based on translation rate and available mRNA
    double fProteinAmount = m_fTranslationRate * dt * m_fNumber;

    // Get sequence from GeneWiki
    const std::string& sequence = GeneWiki::getInstance().getSequence(m_sGeneName);

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
    auto pProtein = std::make_shared<MPopulation>(m_sProteinName, fProteinAmount);

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

void MRNA::splice()
{
    // TODO: Implement RNA splicing
    // This would remove introns and join exons
    // For now, we assume all sequences are pre-spliced
}
