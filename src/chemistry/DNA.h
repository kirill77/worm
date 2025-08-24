#pragma once

#include <vector>
#include <memory>
#include <string>
#include <map>
#include "Gene.h"
#include "MRNA.h"
#include "StringDict.h"

// Forward declarations  
class GridCell;

class DNA
{
private:
    std::vector<std::shared_ptr<Gene>> m_pGenes;
    std::map<StringDict::ID, std::shared_ptr<Gene>> m_geneMap; // Quick lookup by ID

public:
    // Add a gene to the DNA
    void addGene(StringDict::ID id, double expressionRate = 1.0, double basalLevel = 0.1);

    // Get a gene by ID
    std::shared_ptr<Gene> getGene(StringDict::ID id) const;

    // Transcribe all genes
    std::vector<std::shared_ptr<MRNA>> transcribeAll(double dt) const;

    // Regulate gene expression
    void regulateGene(StringDict::ID id, double newExpressionRate);
    
    // Update gene expression based on transcription factors (protein concentrations)
    void updateTranscriptionalRegulation(double dt, const class GridCell& nuclearCompartment);
};

