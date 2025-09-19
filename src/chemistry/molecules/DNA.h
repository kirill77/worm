#pragma once

#include <vector>
#include <memory>
#include <string>
#include <map>
#include "Gene.h"

#include "StringDict.h"

// Forward declarations  
class GridCell;

class DNA
{
private:
    std::vector<std::shared_ptr<Gene>> m_pGenes;
    std::map<StringDict::ID, std::shared_ptr<Gene>> m_geneMap; // Quick lookup by ID
    Species m_species = Species::GENERIC;

public:
    DNA() = default;
    explicit DNA(Species species) : m_species(species) {}

    // Add a gene to the DNA
    void addGene(StringDict::ID id, double expressionRate = 1.0, double basalLevel = 0.1);

    Species getSpecies() const { return m_species; }

    // Get a gene by ID
    std::shared_ptr<Gene> getGene(StringDict::ID id) const;

    // Transcribe all genes
    std::vector<std::shared_ptr<MPopulation>> transcribeAll(double dt) const;

    // Regulate gene expression
    void regulateGene(StringDict::ID id, double newExpressionRate);
    
    // Update gene expression based on transcription factors (protein concentrations)
    void updateTranscriptionalRegulation(double dt, const class GridCell& nuclearCompartment);
};

