#include "pch.h"
#include "ProteinInteractionLoader.h"
#include "StringDict.h"
#include "log/ILog.h"
#include "fileUtils/fileUtils.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>

// Define the expected CSV file names
static const std::string PHOSPHORYLATION_FILE = "phosphorylation.csv";
static const std::string DEPHOSPHORYLATION_FILE = "dephosphorylation.csv";
static const std::string COMPLEX_FORMATION_FILE = "complex_formation.csv";

// Helper function to trim whitespace from a string
static std::string TrimWhitespace(const std::string& str) {
    auto start = std::find_if_not(str.begin(), str.end(), [](unsigned char c) {
        return std::isspace(c);
    });
    
    auto end = std::find_if_not(str.rbegin(), str.rend(), [](unsigned char c) {
        return std::isspace(c);
    }).base();
    
    return (start < end) ? std::string(start, end) : std::string();
}

// Helper function to validate that a protein name exists in StringDict
// Handles complex proteins separated by colons (e.g., "PAR-1:CDC-42")
static void ValidateProteinName(const std::string& proteinName, const std::string& context) {
    if (proteinName.empty()) {
        LOG_ERROR("Empty protein name found in %s", context.c_str());
        assert(false && "Empty protein name in CSV file");
        return;
    }
    
    // Split protein name by colons to handle complex proteins
    std::stringstream ss(proteinName);
    std::string individualProtein;
    
    while (std::getline(ss, individualProtein, ':')) {
        // Trim whitespace from individual protein name
        individualProtein = TrimWhitespace(individualProtein);
        
        if (individualProtein.empty()) {
            LOG_ERROR("Empty individual protein name in complex protein '%s' (context: %s)", 
                      proteinName.c_str(), context.c_str());
            assert(false && "Empty individual protein name in complex protein from CSV file");
            continue;
        }
        
        StringDict::ID id = StringDict::stringToId(individualProtein);
        if (id == StringDict::ID::eUNKNOWN) {
            LOG_ERROR("Individual protein '%s' from complex protein '%s' not found in StringDict (context: %s). This indicates a typo or missing definition in StringDict.", 
                      individualProtein.c_str(), proteinName.c_str(), context.c_str());
            assert(false && "Individual protein name from CSV file not found in StringDict");
        }
    }
}

std::vector<std::shared_ptr<ProteinInteraction>> ProteinInteractionLoader::LoadAllInteractions(const std::string& basePath)
{
    std::vector<std::shared_ptr<ProteinInteraction>> allInteractions;
    
    // Verify that the base path exists
    std::filesystem::path path(basePath);
    if (!std::filesystem::exists(path)) {
        LOG_ERROR("Interaction data directory not found: %s", basePath.c_str());
        return allInteractions;
    }
    
    // Make sure path ends with separator
    std::string directory = basePath;
    if (!directory.empty() && directory.back() != '/' && directory.back() != '\\') {
        directory += '/';
    }
    
    // Load phosphorylation interactions
    std::string phosphorylationPath = directory + PHOSPHORYLATION_FILE;
    if (FileExists(phosphorylationPath)) {
        auto phosphoInteractions = LoadPhosphorylationInteractions(phosphorylationPath);
        LOG_INFO("Loaded %zu phosphorylation interactions from %s", phosphoInteractions.size(), phosphorylationPath.c_str());
        
        // Add to master list
        allInteractions.insert(allInteractions.end(), 
                              phosphoInteractions.begin(), 
                              phosphoInteractions.end());
    }
    
    // Load dephosphorylation interactions
    std::string dephosphorylationPath = directory + DEPHOSPHORYLATION_FILE;
    if (FileExists(dephosphorylationPath)) {
        auto dephosphoInteractions = LoadDephosphorylationInteractions(dephosphorylationPath);
        LOG_INFO("Loaded %zu dephosphorylation interactions from %s", dephosphoInteractions.size(), dephosphorylationPath.c_str());
        
        // Add to master list
        allInteractions.insert(allInteractions.end(), 
                              dephosphoInteractions.begin(), 
                              dephosphoInteractions.end());
    }
    
    // Load complex formation interactions
    std::string complexFormationPath = directory + COMPLEX_FORMATION_FILE;
    if (FileExists(complexFormationPath)) {
        auto complexInteractions = LoadComplexFormationInteractions(complexFormationPath);
        LOG_INFO("Loaded %zu complex formation interactions from %s", complexInteractions.size(), complexFormationPath.c_str());
        
        // Add to master list
        allInteractions.insert(allInteractions.end(), 
                              complexInteractions.begin(), 
                              complexInteractions.end());
    }
    
    return allInteractions;
}

std::vector<std::shared_ptr<PhosphorylationInteraction>> ProteinInteractionLoader::LoadPhosphorylationInteractions(const std::string& filePath)
{
    std::vector<std::shared_ptr<PhosphorylationInteraction>> interactions;
    
    // Open CSV file
    std::ifstream file(filePath);
    if (!file.is_open()) {
        LOG_ERROR("Failed to open phosphorylation interaction file: %s", filePath.c_str());
        return interactions;
    }
    
    // Read header line
    std::string header;
    std::getline(file, header);
    
    // Read and parse each line
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines and comments
        }
        
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> values;
        
        // Parse CSV row
        while (std::getline(ss, cell, ',')) {
            values.push_back(TrimWhitespace(cell));
        }
        
        // Ensure we have enough values
        if (values.size() < 4) {
            LOG_WARN("Skipping malformed phosphorylation entry: %s", line.c_str());
            continue;
        }
        
        try {
            // Extract values
            std::string kinaseName = values[0]; // Already trimmed
            std::string targetName = values[1]; // Already trimmed
            double removalRate = std::stod(values[2]);
            double saturationConstant = std::stod(values[3]);
            
            // Validate protein names
            ValidateProteinName(kinaseName, "phosphorylation kinase");
            ValidateProteinName(targetName, "phosphorylation target");
            
            // Create interaction parameters
            PhosphorylationInteraction::Parameters params{
                removalRate,
                saturationConstant
            };
            
            // Create and add interaction
            interactions.push_back(std::make_shared<PhosphorylationInteraction>(
                kinaseName, targetName, params));
        }
        catch (const std::exception& e) {
            LOG_ERROR("Error parsing phosphorylation interaction: %s - %s", line.c_str(), e.what());
        }
    }
    
    return interactions;
}

std::vector<std::shared_ptr<DephosphorylationInteraction>> ProteinInteractionLoader::LoadDephosphorylationInteractions(const std::string& filePath)
{
    std::vector<std::shared_ptr<DephosphorylationInteraction>> interactions;
    
    // Open CSV file
    std::ifstream file(filePath);
    if (!file.is_open()) {
        LOG_ERROR("Failed to open dephosphorylation interaction file: %s", filePath.c_str());
        return interactions;
    }
    
    // Read header line
    std::string header;
    std::getline(file, header);
    
    // Read and parse each line
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines and comments
        }
        
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> values;
        
        // Parse CSV row
        while (std::getline(ss, cell, ',')) {
            values.push_back(TrimWhitespace(cell));
        }
        
        // Ensure we have enough values
        if (values.size() < 2) {
            LOG_WARN("Skipping malformed dephosphorylation entry: %s", line.c_str());
            continue;
        }
        
        try {
            // Extract values
            std::string targetName = values[0]; // Already trimmed
            double recoveryRate = std::stod(values[1]);
            
            // Validate protein name
            ValidateProteinName(targetName, "dephosphorylation target");
            
            // Create interaction parameters
            DephosphorylationInteraction::Parameters params{
                recoveryRate
            };
            
            // Create and add interaction
            interactions.push_back(std::make_shared<DephosphorylationInteraction>(
                targetName, params));
        }
        catch (const std::exception& e) {
            LOG_ERROR("Error parsing dephosphorylation interaction: %s - %s", line.c_str(), e.what());
        }
    }
    
    return interactions;
}

std::vector<std::shared_ptr<ComplexFormationInteraction>> ProteinInteractionLoader::LoadComplexFormationInteractions(const std::string& filePath)
{
    std::vector<std::shared_ptr<ComplexFormationInteraction>> interactions;
    
    // Open CSV file
    std::ifstream file(filePath);
    if (!file.is_open()) {
        LOG_ERROR("Failed to open complex formation interaction file: %s", filePath.c_str());
        return interactions;
    }
    
    // Read header line
    std::string header;
    std::getline(file, header);
    
    // Read and parse each line
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines and comments
        }
        
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> values;
        
        // Parse CSV row
        while (std::getline(ss, cell, ',')) {
            values.push_back(TrimWhitespace(cell));
        }
        
        // Ensure we have enough values
        if (values.size() < 5) {
            LOG_WARN("Skipping malformed complex formation entry: %s", line.c_str());
            continue;
        }
        
        try {
            // Extract values
            std::string firstProtein = values[0]; // Already trimmed
            std::string secondProtein = values[1]; // Already trimmed
            double bindingRate = std::stod(values[2]);
            double dissociationRate = std::stod(values[3]);
            double saturationConstant = std::stod(values[4]);
            std::string complexName = values.size() > 5 ? values[5] : firstProtein + "-" + secondProtein;
            
            // Validate protein names
            ValidateProteinName(firstProtein, "complex formation first protein");
            ValidateProteinName(secondProtein, "complex formation second protein");
            ValidateProteinName(complexName, "complex formation complex name");
            
            // Create interaction parameters
            ComplexFormationInteraction::Parameters params{
                bindingRate,
                dissociationRate,
                saturationConstant,
                complexName
            };
            
            // Create and add interaction
            interactions.push_back(std::make_shared<ComplexFormationInteraction>(
                firstProtein, secondProtein, params));
        }
        catch (const std::exception& e) {
            LOG_ERROR("Error parsing complex formation interaction: %s - %s", line.c_str(), e.what());
        }
    }
    
    return interactions;
}

bool ProteinInteractionLoader::FileExists(const std::string& filePath)
{
    return std::filesystem::exists(filePath);
} 