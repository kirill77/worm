#pragma once

#include <string>
#include <map>
#include <memory>
#include <vector>
#include <filesystem>
#include "chemistry/molecules/Molecule.h"

class GeneWiki
{
private:
    mutable std::map<std::string, std::string> m_sequences;  // Map of gene name to sequence (mutable for lazy-loading)
    struct GeneData
    {
        std::vector<std::pair<Molecule, uint32_t>> m_trnaRequirements;  // charged tRNA requirements per protein
    };
    mutable std::map<std::string, GeneData> m_geneData;  // Map of gene name to precomputed tRNA requirements

    // IO helpers
    std::filesystem::path getGenesFolder() const;
    std::filesystem::path getGeneFilePath(const std::string& geneName) const;
    static std::string sanitizeGeneNameForFile(const std::string& geneName);
    bool loadSequenceFromFile(const std::filesystem::path& filePath, std::string& outSequence) const;
    bool saveSequenceToFile(const std::filesystem::path& filePath, const std::string& sequence) const;
    bool fetchSequenceFromPublicDb(const std::string& geneName, std::string& outSequence) const;
    bool ensureSequenceLoaded(const std::string& geneName) const;
    bool ensureGeneDataComputed(const std::string& geneName) const;
    void initializeDefaultGeneData();
    
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

    // Get precomputed GeneData for a gene (computed on-demand from sequence)
    const std::vector<std::pair<Molecule, uint32_t>>& getGeneData(const std::string& geneName) const;
    
    // Check if GeneData exists for a gene
    bool hasGeneData(const std::string& geneName) const;
    
}; 