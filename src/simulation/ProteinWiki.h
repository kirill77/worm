#pragma once

#include <vector>
#include <string>
#include <memory>
#include "ProteinAntagonism.h"

// A static repository of protein interaction data
class ProteinWiki
{
private:
    // List of known protein antagonistic relationships
    static std::vector<ProteinAntagonism> s_proteinAntagonisms;

    // Private constructor to prevent instantiation
    ProteinWiki() = default;

public:
    // Initialize all known protein antagonisms
    static void Initialize();

    // Get all known protein antagonisms
    static const std::vector<ProteinAntagonism>& GetProteinAntagonisms();

    // Get protein antagonisms where the specified protein is the antagonist
    static std::vector<ProteinAntagonism> GetAntagonismsFromProtein(const std::string& antagonistName);
    
    // Get protein antagonisms where the specified protein is the target
    static std::vector<ProteinAntagonism> GetAntagonismsToProtein(const std::string& targetName);
}; 