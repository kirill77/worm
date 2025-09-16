#pragma once

#include <string>
#include <map>
#include <memory>
#include <vector>
#include "chemistry/molecules/Molecule.h"

class GeneWiki
{
private:
    std::map<std::string, std::string> m_sequences;  // Map of gene name to sequence
    struct GeneData
    {
        std::vector<std::pair<Molecule, uint32_t>> m_trnaRequirements;  // charged tRNA requirements per protein
    };
    std::map<std::string, GeneData> m_geneData;  // Map of gene name to precomputed tRNA requirements
    
    // Private constructor for singleton pattern
    GeneWiki();

public:
    // Singleton access
    static GeneWiki& getInstance();

    // Delete copy constructor and assignment operator
    GeneWiki(const GeneWiki&) = delete;
    GeneWiki& operator=(const GeneWiki&) = delete;

    // Get sequence for a gene
    const std::string& getSequence(const std::string& geneName) const;
    
    // Check if a sequence exists for a gene
    bool hasSequence(const std::string& geneName) const;

    // Get precomputed GeneData for a gene
    const std::vector<std::pair<Molecule, uint32_t>>& getGeneData(const std::string& geneName) const;
    
    // Check if GeneData exists for a gene
    bool hasGeneData(const std::string& geneName) const;

    // Initialize with default sequences
    void initializeDefaultSequences();
    void initializeDefaultGeneData();
}; 