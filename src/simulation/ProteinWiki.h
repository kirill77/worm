#pragma once

#include <vector>
#include <string>
#include <memory>
#include "ProteinInteraction.h"

// A static repository of protein interaction data
class ProteinWiki
{
private:
    // List of known protein interactions
    static std::vector<std::shared_ptr<ProteinInteraction>> s_proteinInteractions;

    // Private constructor to prevent instantiation
    ProteinWiki() = default;

public:
    // Initialize all known protein interactions
    static void Initialize();

    // Get all known protein interactions
    static const std::vector<std::shared_ptr<ProteinInteraction>>& GetProteinInteractions();

    // Get interactions by mechanism
    static std::vector<std::shared_ptr<ProteinInteraction>> GetInteractionsByMechanism(ProteinInteraction::Mechanism mechanism);
    
    // Utility function to get the name of a phosphorylated protein
    static std::string GetPhosphorylatedName(const std::string& proteinName);
}; 