#include "pch.h"
#include "DNA.h"

void DNA::addGene(const std::string& name, double expressionRate, double basalLevel)
{
    auto gene = std::make_shared<Gene>(name, expressionRate, basalLevel);
    m_pGenes.push_back(gene);
    m_geneMap[name] = gene;
}

std::shared_ptr<Gene> DNA::getGene(const std::string& name) const
{
    auto it = m_geneMap.find(name);
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

void DNA::regulateGene(const std::string& name, double newExpressionRate)
{
    if (auto gene = getGene(name))
    {
        gene->setExpressionRate(newExpressionRate);
    }
}
