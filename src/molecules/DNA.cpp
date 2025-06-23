#include "pch.h"
#include "DNA.h"
#include "StringDict.h"

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
