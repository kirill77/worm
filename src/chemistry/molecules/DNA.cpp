#include "DNA.h"
#include "StringDict.h"
#include "GridCell.h"
#include "TRNA.h"
#include "simConstants.h"

void DNA::addGene(StringDict::ID id, double expressionRate, double basalLevel)
{
    auto gene = std::make_shared<Gene>(id, expressionRate, basalLevel, m_species);
    m_pGenes.push_back(gene);
    m_geneMap[id] = gene;
}

std::shared_ptr<Gene> DNA::getGene(StringDict::ID id) const
{
    auto it = m_geneMap.find(id);
    return (it != m_geneMap.end()) ? it->second : nullptr;
}

std::vector<std::shared_ptr<MPopulation>> DNA::transcribeAll(double dt) const
{
    std::vector<std::shared_ptr<MPopulation>> transcribedRNA;
    
    for (const auto& gene : m_pGenes)
    {
        // For tRNA genes, produce TRNA molecules directly (Pol III products), not mRNA
        StringDict::ID id = gene->getId();
        bool isTRNAGene = TRNA::isTRNAGeneId(id);

        if (isTRNAGene)
        {
            if (auto amountProbe = gene->transcribe(dt))
            {
                // Reuse computed production amount; emit uncharged TRNA of the same gene ID and species
                // TRNAs are species-agnostic in this model; use GENERIC to match charging/translation keys
                Molecule trnaMol(id, ChemicalType::TRNA, Species::GENERIC);
                // Diagnostic boost: increase tRNA nuclear production to test charging/consumption bottlenecks
                const double kTRNA_PRODUCTION_BOOST = MoleculeConstants::TRNA_POLIII_PRODUCTION_MULTIPLIER;
                Population boosted = amountProbe->m_population;
                boosted.m_fNumber *= kTRNA_PRODUCTION_BOOST;
                auto trnaPop = std::make_shared<MPopulation>(trnaMol, boosted);
                transcribedRNA.push_back(trnaPop);
            }
        }
        else
        {
            if (auto mRNA = gene->transcribe(dt))
            {
                transcribedRNA.push_back(mRNA);
            }
        }
    }
    
    return transcribedRNA;
}

void DNA::regulateGene(StringDict::ID id, double newExpressionRate)
{
    if (auto gene = getGene(id))
    {
        gene->setExpressionRate(newExpressionRate);
    }
}

void DNA::updateTranscriptionalRegulation(double dt, const GridCell& nuclearCompartment)
{
    // Regulate γ-tubulin gene expression based on CDK2/CyclinE levels
    // This mimics E2F transcription factor activity during S/G2 phases
    if (auto gammaTubulinGene = getGene(StringDict::ID::GAMMA_TUBULIN))
    {
        // Only do expensive protein lookups if gene exists
        double cdk2Level = 0.0;
        auto cdk2It = nuclearCompartment.m_molecules.find(Molecule(StringDict::ID::CDK_2, ChemicalType::PROTEIN));
        if (cdk2It != nuclearCompartment.m_molecules.end()) {
            cdk2Level = cdk2It->second.m_fNumber;
        }
        
        double cyclinELevel = 0.0;
        auto cyclinEIt = nuclearCompartment.m_molecules.find(Molecule(StringDict::ID::CCE_1, ChemicalType::PROTEIN));
        if (cyclinEIt != nuclearCompartment.m_molecules.end()) {
            cyclinELevel = cyclinEIt->second.m_fNumber;
        }
        
        // Calculate transcriptional activation using Hill kinetics
        // Both CDK2 and CyclinE needed for activation (AND logic)
        double transcriptionFactorActivity = (cdk2Level * cyclinELevel) / 
                                           (MoleculeConstants::TF_ACTIVITY_K + (cdk2Level * cyclinELevel));
        
        // Base expression rate + cell cycle-activated rate
        double basalRate = MoleculeConstants::TRANSCRIPTION_BASAL_RATE;  // Low basal transcription
        double maxActivatedRate = MoleculeConstants::TRANSCRIPTION_MAX_ACTIVATED_RATE;  // Maximum activated transcription
        double newExpressionRate = basalRate + (maxActivatedRate * transcriptionFactorActivity);
        
        // Set the new expression rate
        gammaTubulinGene->setExpressionRate(newExpressionRate);
    }
}