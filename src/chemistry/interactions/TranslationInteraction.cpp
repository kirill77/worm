#include "TranslationInteraction.h"
#include "ResourceDistributor.h"
#include "GridCell.h"
#include "chemistry/molecules/GeneWiki.h"
#include "chemistry/molecules/TRNA.h"
#include "chemistry/molecules/StringDict.h"
#include <algorithm>
#include <cmath>
#include <cassert>

TranslationInteraction::TranslationInteraction(const Molecule& mRNA, const Parameters& params)
    : MoleculeInteraction(Mechanism::TRANSLATION, 0.3)  // ATP cost for translation
    , m_mRNA(mRNA)
    , m_translationRate(params.translationRate)
{
    // Ensure we're dealing with an mRNA molecule
    assert(mRNA.getType() == ChemicalType::MRNA && "TranslationInteraction requires an mRNA molecule");
}

bool TranslationInteraction::apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const
{
    // Get mRNA amount available
    double mRNAAmount = resDistributor.getAvailableResource(m_mRNA);
    
    // Check if we have enough mRNA to produce protein
    if (mRNAAmount < 0.01) {
        return false;  // Not enough mRNA for translation
    }
    
    // Calculate potential protein production
    double potentialProteinAmount = m_translationRate * dt * mRNAAmount;
    
    // Query precomputed tRNA requirements for this gene
    const auto& geneTRNAs = GeneWiki::getInstance().getGeneData(m_mRNA.getName());
    
    // Calculate actual protein amount we can produce based on available resources
    double actualProteinAmount = potentialProteinAmount;
    
    // Check resource availability for tRNAs
    for (const auto& tRNAReq : geneTRNAs) {
        const Molecule& trnaMol = tRNAReq.first;
        uint32_t count = tRNAReq.second;
        if (count == 0) continue;
        double availableTRNA = resDistributor.getAvailableResource(trnaMol);
        double requiredTRNA = static_cast<double>(count) * potentialProteinAmount;
        
        if (availableTRNA < requiredTRNA) {
            // Limit protein production by available tRNA
            actualProteinAmount = std::min(actualProteinAmount, availableTRNA / static_cast<double>(count));
        }
    }
    
    // Calculate the ATP cost
    double requiredATP = actualProteinAmount * m_atpCost;
    
    // If we're in a dry run, just report resource requirements
    if (resDistributor.isDryRun()) {
        if (actualProteinAmount > 0) {
            // Register resource requirements
            resDistributor.notifyResourceWanted(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE), requiredATP);
            resDistributor.notifyResourceWanted(m_mRNA, actualProteinAmount / m_translationRate / dt);
            
            for (const auto& tRNAReq : geneTRNAs) {
                if (tRNAReq.second == 0) continue;
                double requiredTRNA = static_cast<double>(tRNAReq.second) * actualProteinAmount;
                resDistributor.notifyResourceWanted(tRNAReq.first, requiredTRNA);
            }
            return true;
        }
        return false;
    }
    
    // Consume ATP directly from the cell
    auto& atpPop = cell.getOrCreateMolPop(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE));
    if (atpPop.m_fNumber < requiredATP) {
        return false;  // Not enough ATP
    }
    atpPop.m_fNumber -= requiredATP;
    assert(atpPop.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL);
    
    // Don't consume mRNA (it can be translated multiple times)
    // But we do consume tRNAs
    consumeTRNAs(cell, geneTRNAs, actualProteinAmount);
    
    // Create the protein
    Molecule protein(m_mRNA.getID(), ChemicalType::PROTEIN);  // Same ID but different type
    Population& proteinPop = cell.getOrCreateMolPop(protein);
    proteinPop.m_fNumber += actualProteinAmount;
    
    return actualProteinAmount > 0;
}

void TranslationInteraction::consumeTRNAs(GridCell& cell,
                                        const std::vector<std::pair<Molecule, uint32_t>>& requiredTRNAs,
                                        double proteinAmount) const
{
    for (const auto& tRNAReq : requiredTRNAs) {
        auto it = cell.m_molecules.find(tRNAReq.first);
        if (it != cell.m_molecules.end()) {
            double consumeAmount = static_cast<double>(tRNAReq.second) * proteinAmount;
            it->second.m_fNumber -= consumeAmount;
            if (it->second.m_fNumber < 0) {
                it->second.m_fNumber = 0;
            }
        }
    }
}
