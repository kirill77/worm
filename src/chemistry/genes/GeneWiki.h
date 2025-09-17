#pragma once

#include <string>
#include <map>
#include <memory>
#include <vector>
#include <filesystem>
#include <unordered_set>
#include <unordered_map>
#include "chemistry/molecules/Molecule.h"

class GeneWiki
{
private:
    // Keyed by species+gene composite key
    mutable std::map<std::string, std::string> m_sequences;  // Map of (species::gene) to sequence (mutable for lazy-loading)
    struct GeneData
    {
        std::vector<std::pair<Molecule, uint32_t>> m_trnaRequirements;  // charged tRNA requirements per protein
    };
    mutable std::unordered_map<Molecule, GeneData> m_geneData;  // keyed by mRNA molecule (species-aware)
    // Persistent negative cache of missing sequences (keyed by mRNA molecule, species-aware)
    mutable std::unordered_set<Molecule> m_missingSequenceKeys;

    // IO helpers
    std::filesystem::path getGenesFolder() const;
    std::filesystem::path getSpeciesFolder(Species species) const;
    std::filesystem::path getGeneFilePath(Species species, const std::string& geneName) const;
    std::filesystem::path getMissingCacheFilePath(Species species) const;
    static std::string sanitizeGeneNameForFile(const std::string& geneName);
    bool loadSequenceFromFile(const std::filesystem::path& filePath, std::string& outSequence) const;
    bool saveSequenceToFile(const std::filesystem::path& filePath, const std::string& sequence) const;
    bool fetchSequenceFromPublicDb(Species species, const std::string& geneName, std::string& outSequence) const;
    bool ensureSequenceLoaded(Species species, const std::string& geneName) const;
    void loadMissingCache(Species species) const;
    void saveMissingCache(Species species) const;
    bool isMarkedMissing(Species species, const std::string& geneName) const;
    void markMissing(Species species, const std::string& geneName) const;
    void markFound(Species species, const std::string& geneName) const;
    bool ensureGeneDataComputed(Species species, const std::string& geneName) const;
    void initializeDefaultGeneData();
    static std::string makeKey(Species species, const std::string& geneName);
    
    // Private constructor for singleton pattern
    GeneWiki();

public:
    // Singleton access
    static GeneWiki& getInstance();

    // Delete copy constructor and assignment operator
    GeneWiki(const GeneWiki&) = delete;
    GeneWiki& operator=(const GeneWiki&) = delete;

    // Get precomputed GeneData for a gene (computed on-demand from sequence)
    const std::vector<std::pair<Molecule, uint32_t>>& getGeneData(const Molecule& geneMolecule) const;
    
    // Check if GeneData exists for a gene
    bool hasGeneData(const Molecule& geneMolecule) const;
    
}; 