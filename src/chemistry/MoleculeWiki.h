#pragma once

#include <vector>
#include <string>
#include <memory>
#include "ProteinInteraction.h"
#include "StringDict.h"

// A static repository of molecule interaction data
class MoleculeWiki
{
public:
    // Private constructor to prevent instantiation
    MoleculeWiki() = default;
    
private:
    // List of known protein interactions (molecules include proteins)
    static std::vector<std::shared_ptr<ProteinInteraction>> s_proteinInteractions;

    // Helper method to load default hardcoded interactions
    static void LoadDefaultInteractions();

public:
    // Initialize all known molecule interactions
    static void Initialize();

    // Get all known protein interactions
    static const std::vector<std::shared_ptr<ProteinInteraction>>& GetProteinInteractions();

    // Get interactions by mechanism
    static std::vector<std::shared_ptr<ProteinInteraction>> GetInteractionsByMechanism(ProteinInteraction::Mechanism mechanism);
    
    // Utility function to get the name of a phosphorylated protein/molecule
    static std::string GetPhosphorylatedName(const std::string& proteinName);
    
    // Get the bound molecule name for a molecule on a specific surface
    static std::string GetBoundProteinName(const std::string& proteinName, StringDict::ID surface);
};
