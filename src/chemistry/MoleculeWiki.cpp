#include "pch.h"
#include "MoleculeWiki.h"
#include "PhosphorylationInteraction.h"
#include "ComplexFormationInteraction.h"
#include "DephosphorylationInteraction.h"
#include "BindingSurface.h"
#include "ProteinInteractionLoader.h"
#include "utils/log/ILog.h"
#include "utils/fileUtils/fileUtils.h"
#include <algorithm>
#include <iterator>
#include <filesystem>

// Initialize static members
std::vector<std::shared_ptr<ProteinInteraction>> MoleculeWiki::s_proteinInteractions;
std::unordered_map<Molecule, MolInfo> MoleculeWiki::m_moleculesInfo;

void MoleculeWiki::Initialize()
{
    // Clear any existing interactions
    s_proteinInteractions.clear();
    
    // Try to find the data directory
    std::filesystem::path dataPath;
    bool dataFolderFound = false;
    
    // First try to find a "data" folder relative to the current directory
    if (std::filesystem::exists("data/proteinRules")) {
        dataPath = "data/proteinRules";
        dataFolderFound = true;
    } 
    // Then try to use the fileUtils helper
    else if (FileUtils::findTheFolder("data", dataPath)) {
        dataPath /= "proteinRules";
        if (std::filesystem::exists(dataPath)) {
            dataFolderFound = true;
        }
    }
    
    // Try a few common paths relative to executable
    if (!dataFolderFound) {
        std::vector<std::string> commonPaths = {
            "../data/proteinRules",
            "../../data/proteinRules",
            "../../../data/proteinRules"
        };
        
        for (const auto& path : commonPaths) {
            if (std::filesystem::exists(path)) {
                dataPath = path;
                dataFolderFound = true;
                break;
            }
        }
    }
    
    if (dataFolderFound) {
        LOG_INFO("Loading molecule interactions from %s", dataPath.string().c_str());
        s_proteinInteractions = ProteinInteractionLoader::LoadAllInteractions(dataPath.string());
        if (s_proteinInteractions.empty()) {
            LOG_ERROR("No molecule interactions were loaded from CSV files.");
        }
    } else {
        LOG_ERROR("Interaction data directory not found. Using default hardcoded interactions.");
    }
}

const std::vector<std::shared_ptr<ProteinInteraction>>& MoleculeWiki::GetProteinInteractions()
{
    return s_proteinInteractions;
}

std::vector<std::shared_ptr<ProteinInteraction>> MoleculeWiki::GetInteractionsByMechanism(ProteinInteraction::Mechanism mechanism)
{
    std::vector<std::shared_ptr<ProteinInteraction>> result;
    
    std::copy_if(s_proteinInteractions.begin(), s_proteinInteractions.end(), 
                 std::back_inserter(result),
                 [mechanism](const auto& interaction) {
                     return interaction->getMechanism() == mechanism;
                 });
                 
    return result;
}

std::string MoleculeWiki::GetPhosphorylatedName(const std::string& proteinName)
{
    return proteinName + "-P";
}

std::string MoleculeWiki::GetBoundProteinName(const std::string& proteinName, StringDict::ID surface)
{
    return proteinName + ":" + StringDict::idToString(surface);
}

const MolInfo& MoleculeWiki::getInfo(const Molecule& molecule)
{
    auto it = m_moleculesInfo.find(molecule);
    if (it != m_moleculesInfo.end()) {
        return it->second;
    }
    
    // If not found, return a default empty MolInfo object
    static const MolInfo defaultInfo("No information available", "", 0.0, "", 0.0, 0.0);
    return defaultInfo;
}
