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
    
    // Get sequence from GeneWiki
    const std::string& sequence = GeneWiki::getInstance().getSequence(m_mRNA.getName());
    
    // Check for required tRNAs and calculate what we actually need
    std::vector<std::pair<Molecule, double>> requiredTRNAs;
    if (!checkRequiredTRNAs(cell.m_molecules, sequence, requiredTRNAs)) {
        return false;  // Missing required tRNAs
    }
    
    // Calculate actual protein amount we can produce based on available resources
    double actualProteinAmount = potentialProteinAmount;
    
    // Check resource availability for tRNAs
    for (const auto& tRNAReq : requiredTRNAs) {
        double availableTRNA = resDistributor.getAvailableResource(tRNAReq.first);
        double requiredTRNA = tRNAReq.second * potentialProteinAmount;
        
        if (availableTRNA < requiredTRNA) {
            // Limit protein production by available tRNA
            actualProteinAmount = std::min(actualProteinAmount, availableTRNA / tRNAReq.second);
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
            
            for (const auto& tRNAReq : requiredTRNAs) {
                double requiredTRNA = tRNAReq.second * actualProteinAmount;
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
    consumeTRNAs(cell, requiredTRNAs, actualProteinAmount);
    
    // Create the protein
    Molecule protein(m_mRNA.getID(), ChemicalType::PROTEIN);  // Same ID but different type
    Population& proteinPop = cell.getOrCreateMolPop(protein);
    proteinPop.m_fNumber += actualProteinAmount;
    
    return actualProteinAmount > 0;
}

bool TranslationInteraction::checkRequiredTRNAs(const std::unordered_map<Molecule, Population>& availableMolecules,
                                               const std::string& sequence,
                                               std::vector<std::pair<Molecule, double>>& requiredTRNAs) const
{
    requiredTRNAs.clear();
    
    // Track unique tRNAs needed
    std::unordered_map<StringDict::ID, double> tRNAUsage;
    
    // Check each codon in the sequence
    for (size_t i = 0; i < sequence.length(); i += 3) {
        if (i + 2 >= sequence.length()) break;  // Incomplete codon
        
        std::string codon = sequence.substr(i, 3);
        bool codonMatched = false;
        
        // Convert codon to anticodon and find matching charged tRNAs
        std::string anticodon = TRNA::codonToAnticodon(codon);
        if (!anticodon.empty()) {
            auto matchingTRNAIds = TRNA::getChargedTRNAsWithAnticodon(anticodon);
            
            for (StringDict::ID tRNAId : matchingTRNAIds) {
                Molecule tRNAMolecule(tRNAId, ChemicalType::TRNA);
                auto it = availableMolecules.find(tRNAMolecule);
                if (it != availableMolecules.end() && it->second.m_fNumber > 0.01) {
                    codonMatched = true;
                    tRNAUsage[tRNAId] += 0.01;  // Each codon consumes 0.01 units of tRNA
                    break;
                }
            }
        }
        
        if (!codonMatched) {
            return false;  // Missing required tRNA
        }
    }
    
    // Convert usage map to required tRNAs vector
    for (const auto& usage : tRNAUsage) {
        Molecule tRNA(usage.first, ChemicalType::TRNA);
        requiredTRNAs.emplace_back(tRNA, usage.second);
    }
    
    return true;
}

void TranslationInteraction::consumeTRNAs(GridCell& cell,
                                        const std::vector<std::pair<Molecule, double>>& requiredTRNAs,
                                        double proteinAmount) const
{
    for (const auto& tRNAReq : requiredTRNAs) {
        auto it = cell.m_molecules.find(tRNAReq.first);
        if (it != cell.m_molecules.end()) {
            double consumeAmount = tRNAReq.second * proteinAmount;
            it->second.m_fNumber -= consumeAmount;
            if (it->second.m_fNumber < 0) {
                it->second.m_fNumber = 0;
            }
        }
    }
}
