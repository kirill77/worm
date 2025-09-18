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
    // Species-specific alias for public DB lookup: map mRNA molecule → canonical symbol/ID
    mutable std::unordered_map<Molecule, std::string> m_lookupAliases;

    // IO helpers
    std::filesystem::path getGenesFolder() const;
    std::filesystem::path getSpeciesFolder(Species species) const;
    std::filesystem::path getGeneFilePath(Species species, const std::string& geneName) const;
    std::filesystem::path getMissingCacheFilePath(Species species) const;
    static std::string sanitizeGeneNameForFile(const std::string& geneName);
    bool loadSequenceFromFile(const std::filesystem::path& filePath, std::string& outSequence) const;
    bool saveSequenceToFile(const std::filesystem::path& filePath, const std::string& sequence) const;
    // Public DB fetch now takes Molecule for species-safe lookup
    bool fetchSequenceFromPublicDb(const Molecule& mrna, std::string& outSequence) const;
    // Resolve species-specific lookup name; Molecule overload is the only variant
    std::string resolveLookupName(const Molecule& mrna) const;
    bool ensureSequenceLoaded(const Molecule& mrna) const;
    void loadMissingCache(Species species) const;
    void saveMissingCache(Species species) const;
    bool isMarkedMissing(const Molecule& mrna) const;
    void markMissing(const Molecule& mrna) const;
    void markFound(const Molecule& mrna) const;
    bool ensureGeneDataComputed(const Molecule& mrna) const;

    static std::string makeKey(Species species, const std::string& geneName);

    // Helper: map codon (3 chars) to charged tRNA ID
    static StringDict::ID codonToChargedTrnaId(const std::string &codon);

public:
    GeneWiki();
    static GeneWiki& getInstance();

    // Accessors for gene data (charged tRNA requirements per protein)
    const std::vector<std::pair<Molecule, uint32_t>>& getGeneData(const Molecule& geneMolecule) const;
    bool hasGeneData(const Molecule& geneMolecule) const;
}; 