#pragma once

#include <string>
#include <map>
#include <memory>

class GeneWiki
{
private:
    std::map<std::string, std::string> m_sequences;  // Map of gene name to sequence
    
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

    // Initialize with default sequences
    void initializeDefaultSequences();
}; 