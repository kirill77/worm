#pragma once

#include <vector>
#include <string>
#include <memory>
#include "ProteinInteraction.h"
#include "StringDict.h"

// A static repository of protein interaction data
class ProteinWiki
{
public:
    // Private constructor to prevent instantiation
    ProteinWiki() = default;
    
private:
    // List of known protein interactions
    static std::vector<std::shared_ptr<ProteinInteraction>> s_proteinInteractions;

    // Helper method to load default hardcoded interactions
    static void LoadDefaultInteractions();

public:
    // Initialize all known protein interactions
    static void Initialize();

    // Get all known protein interactions
    static const std::vector<std::shared_ptr<ProteinInteraction>>& GetProteinInteractions();

    // Get interactions by mechanism
    static std::vector<std::shared_ptr<ProteinInteraction>> GetInteractionsByMechanism(ProteinInteraction::Mechanism mechanism);
    
    // Utility function to get the name of a phosphorylated protein
    static std::string GetPhosphorylatedName(const std::string& proteinName);
    
    // Get the binding site name for a specific surface
    static std::string GetBindingSiteName(StringDict::ID surface);
    
    // Get the bound protein name for a protein on a specific surface
    static std::string GetBoundProteinName(const std::string& proteinName, StringDict::ID surface);
}; 