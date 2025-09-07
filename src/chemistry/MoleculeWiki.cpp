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
    m_moleculesInfo.clear();
    
    // Initialize tRNA molecule information
    initializeTRNAInfo();
    
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

void MoleculeWiki::initializeTRNAInfo()
{
    // Initialize uncharged tRNA molecules with their charging rates
    
    // Start codon - absolutely essential
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_MET_ATG, ChemicalType::TRNA)] = 
        MolInfo("Methionine tRNA (uncharged)", "tRNA", 25000.0, "Start codon tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_MET_ATG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Methionine tRNA (charged)", "tRNA", 25000.0, "Start codon tRNA", 0.0, 0.0, 0.0);
    
    // Common amino acids - high charging rates
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLY_GGA, ChemicalType::TRNA)] = 
        MolInfo("Glycine tRNA GGA (uncharged)", "tRNA", 25000.0, "Glycine tRNA", 0.0, 0.0, 0.9);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLY_GGA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Glycine tRNA GGA (charged)", "tRNA", 25000.0, "Glycine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLY_GGT, ChemicalType::TRNA)] = 
        MolInfo("Glycine tRNA GGT (uncharged)", "tRNA", 25000.0, "Glycine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLY_GGT_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Glycine tRNA GGT (charged)", "tRNA", 25000.0, "Glycine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ALA_GCA, ChemicalType::TRNA)] = 
        MolInfo("Alanine tRNA GCA (uncharged)", "tRNA", 25000.0, "Alanine tRNA", 0.0, 0.0, 0.9);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ALA_GCA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Alanine tRNA GCA (charged)", "tRNA", 25000.0, "Alanine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ALA_GCC, ChemicalType::TRNA)] = 
        MolInfo("Alanine tRNA GCC (uncharged)", "tRNA", 25000.0, "Alanine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ALA_GCC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Alanine tRNA GCC (charged)", "tRNA", 25000.0, "Alanine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LEU_CTG, ChemicalType::TRNA)] = 
        MolInfo("Leucine tRNA CTG (uncharged)", "tRNA", 25000.0, "Leucine tRNA", 0.0, 0.0, 0.9);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LEU_CTG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Leucine tRNA CTG (charged)", "tRNA", 25000.0, "Leucine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LEU_CTC, ChemicalType::TRNA)] = 
        MolInfo("Leucine tRNA CTC (uncharged)", "tRNA", 25000.0, "Leucine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LEU_CTC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Leucine tRNA CTC (charged)", "tRNA", 25000.0, "Leucine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_SER_TCA, ChemicalType::TRNA)] = 
        MolInfo("Serine tRNA TCA (uncharged)", "tRNA", 25000.0, "Serine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_SER_TCA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Serine tRNA TCA (charged)", "tRNA", 25000.0, "Serine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_SER_TCG, ChemicalType::TRNA)] = 
        MolInfo("Serine tRNA TCG (uncharged)", "tRNA", 25000.0, "Serine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_SER_TCG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Serine tRNA TCG (charged)", "tRNA", 25000.0, "Serine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_VAL_GTG, ChemicalType::TRNA)] = 
        MolInfo("Valine tRNA GTG (uncharged)", "tRNA", 25000.0, "Valine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_VAL_GTG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Valine tRNA GTG (charged)", "tRNA", 25000.0, "Valine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_VAL_GTC, ChemicalType::TRNA)] = 
        MolInfo("Valine tRNA GTC (uncharged)", "tRNA", 25000.0, "Valine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_VAL_GTC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Valine tRNA GTC (charged)", "tRNA", 25000.0, "Valine tRNA", 0.0, 0.0, 0.0);
        
    // Essential amino acids - lower charging rates
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LYS_AAG, ChemicalType::TRNA)] = 
        MolInfo("Lysine tRNA AAG (uncharged)", "tRNA", 25000.0, "Lysine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LYS_AAG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Lysine tRNA AAG (charged)", "tRNA", 25000.0, "Lysine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ASP_GAC, ChemicalType::TRNA)] = 
        MolInfo("Aspartic acid tRNA GAC (uncharged)", "tRNA", 25000.0, "Aspartic acid tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ASP_GAC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Aspartic acid tRNA GAC (charged)", "tRNA", 25000.0, "Aspartic acid tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLU_GAG, ChemicalType::TRNA)] = 
        MolInfo("Glutamic acid tRNA GAG (uncharged)", "tRNA", 25000.0, "Glutamic acid tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLU_GAG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Glutamic acid tRNA GAG (charged)", "tRNA", 25000.0, "Glutamic acid tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_PRO_CCA, ChemicalType::TRNA)] = 
        MolInfo("Proline tRNA CCA (uncharged)", "tRNA", 25000.0, "Proline tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_PRO_CCA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Proline tRNA CCA (charged)", "tRNA", 25000.0, "Proline tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_THR_ACA, ChemicalType::TRNA)] = 
        MolInfo("Threonine tRNA ACA (uncharged)", "tRNA", 25000.0, "Threonine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_THR_ACA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Threonine tRNA ACA (charged)", "tRNA", 25000.0, "Threonine tRNA", 0.0, 0.0, 0.0);
        
    // Add remaining tRNAs with standard rates...
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ARG_CGA, ChemicalType::TRNA)] = 
        MolInfo("Arginine tRNA CGA (uncharged)", "tRNA", 25000.0, "Arginine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ARG_CGA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Arginine tRNA CGA (charged)", "tRNA", 25000.0, "Arginine tRNA", 0.0, 0.0, 0.0);
}