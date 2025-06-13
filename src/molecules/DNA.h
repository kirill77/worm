#pragma once

#include <vector>
#include <memory>
#include <string>
#include <map>
#include "Gene.h"
#include "MRNA.h"

class DNA
{
private:
    std::vector<std::shared_ptr<Gene>> m_pGenes;
    std::map<std::string, std::shared_ptr<Gene>> m_geneMap; // Quick lookup by name

public:
    // Add a gene to the DNA
    void addGene(const std::string& name, double expressionRate = 1.0, double basalLevel = 0.1);

    // Get a gene by name
    std::shared_ptr<Gene> getGene(const std::string& name) const;

    // Transcribe all genes
    std::vector<std::shared_ptr<MRNA>> transcribeAll(double dt) const;

    // Regulate gene expression
    void regulateGene(const std::string& name, double newExpressionRate);
};

