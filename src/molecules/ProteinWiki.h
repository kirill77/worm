#pragma once

#include <vector>
#include <string>
#include <memory>
#include "ProteinInteraction.h"

// A static repository of protein interaction data
class ProteinWiki
{
public:
    // Enum for binding surfaces
    enum class BindingSurface {
        eUNKNOWN,
        MEMBRANE,
        CORTEX,
        CENTROSOME
    };
    
private:
    // List of known protein interactions
    static std::vector<std::shared_ptr<ProteinInteraction>> s_proteinInteractions;

    // Private constructor to prevent instantiation
    ProteinWiki() = default;
    
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
    static std::string GetBindingSiteName(BindingSurface surface);
    
    // Get the bound protein name for a protein on a specific surface
    static std::string GetBoundProteinName(const std::string& proteinName, BindingSurface surface);
    
    // Convert binding surface enum to string
    static std::string BindingSurfaceToString(BindingSurface surface);
}; 