#include "pch.h"
#include "DNA.h"
#include "StringDict.h"
#include "GridCell.h"

void DNA::addGene(StringDict::ID id, double expressionRate, double basalLevel)
{
    auto gene = std::make_shared<Gene>(id, expressionRate, basalLevel);
    m_pGenes.push_back(gene);
    m_geneMap[id] = gene;
}

std::shared_ptr<Gene> DNA::getGene(StringDict::ID id) const
{
    auto it = m_geneMap.find(id);
    return (it != m_geneMap.end()) ? it->second : nullptr;
}

std::vector<std::shared_ptr<MRNA>> DNA::transcribeAll(double dt) const
{
    std::vector<std::shared_ptr<MRNA>> transcribedRNA;
    
    for (const auto& gene : m_pGenes)
    {
        if (auto mRNA = gene->transcribe(dt))
        {
            transcribedRNA.push_back(mRNA);
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
                                           (250000.0 + (cdk2Level * cyclinELevel));
        
        // Base expression rate + cell cycle-activated rate
        double basalRate = 0.05;  // Low basal transcription
        double maxActivatedRate = 0.8;  // Maximum activated transcription
        double newExpressionRate = basalRate + (maxActivatedRate * transcriptionFactorActivity);
        
        // Set the new expression rate
        gammaTubulinGene->setExpressionRate(newExpressionRate);
    }
}