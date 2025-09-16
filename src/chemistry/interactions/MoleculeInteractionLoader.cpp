#include "MoleculeInteractionLoader.h"
#include "chemistry/molecules/StringDict.h"
#include "chemistry/molecules/GeneWiki.h"
#include "chemistry/molecules/MoleculeWiki.h"
#include "utils/log/ILog.h"
#include "utils/fileUtils/fileUtils.h"
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

std::vector<std::shared_ptr<MoleculeInteraction>> MoleculeInteractionLoader::LoadAllInteractions(const std::string& basePath)
{
    std::vector<std::shared_ptr<MoleculeInteraction>> allInteractions;
    
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
    
    // Load translation interactions for all mRNA molecules
    auto translationInteractions = LoadTranslationInteractions();
    LOG_INFO("Created %zu translation interactions", translationInteractions.size());
    
    // Add to master list
    allInteractions.insert(allInteractions.end(),
                          translationInteractions.begin(),
                          translationInteractions.end());
    
    return allInteractions;
}

std::vector<std::shared_ptr<PhosphorylationInteraction>> MoleculeInteractionLoader::LoadPhosphorylationInteractions(const std::string& filePath)
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
        if (values.size() < 5) {
            LOG_WARN("Skipping malformed phosphorylation entry: %s", line.c_str());
            continue;
        }
        
        try {
            // Extract values
            std::string kinaseName = values[0]; // Already trimmed
            std::string targetName = values[1]; // Already trimmed
            std::string phosphorylatedName = values[2]; // Already trimmed
            double removalRate = std::stod(values[3]);
            double saturationConstant = std::stod(values[4]);
            
            // Validate protein names
            ValidateProteinName(kinaseName, "phosphorylation kinase");
            ValidateProteinName(targetName, "phosphorylation target");
            ValidateProteinName(phosphorylatedName, "phosphorylated form");
            
            // Create interaction parameters
            PhosphorylationInteraction::Parameters params{
                removalRate,
                saturationConstant
            };
            
            // Convert string names to StringDict IDs
            StringDict::ID kinaseId = StringDict::stringToId(kinaseName);
            StringDict::ID targetId = StringDict::stringToId(targetName);
            StringDict::ID phosphorylatedId = StringDict::stringToId(phosphorylatedName);
            
            // Create and add interaction
            interactions.push_back(std::make_shared<PhosphorylationInteraction>(
                kinaseId, targetId, phosphorylatedId, params));
        }
        catch (const std::exception& e) {
            LOG_ERROR("Error parsing phosphorylation interaction: %s - %s", line.c_str(), e.what());
        }
    }
    
    return interactions;
}

std::vector<std::shared_ptr<DephosphorylationInteraction>> MoleculeInteractionLoader::LoadDephosphorylationInteractions(const std::string& filePath)
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
        if (values.size() < 3) {
            LOG_WARN("Skipping malformed dephosphorylation entry: %s", line.c_str());
            continue;
        }
        
        try {
            // Extract values
            std::string targetName = values[0]; // Already trimmed
            std::string phosphorylatedName = values[1]; // Already trimmed
            double recoveryRate = std::stod(values[2]);
            
            // Validate protein names
            ValidateProteinName(targetName, "dephosphorylation target");
            ValidateProteinName(phosphorylatedName, "phosphorylated form");
            
            // Create interaction parameters
            DephosphorylationInteraction::Parameters params{
                recoveryRate
            };
            
            // Convert string names to StringDict IDs
            StringDict::ID targetId = StringDict::stringToId(targetName);
            StringDict::ID phosphorylatedId = StringDict::stringToId(phosphorylatedName);
            
            // Create and add interaction
            interactions.push_back(std::make_shared<DephosphorylationInteraction>(
                targetId, phosphorylatedId, params));
        }
        catch (const std::exception& e) {
            LOG_ERROR("Error parsing dephosphorylation interaction: %s - %s", line.c_str(), e.what());
        }
    }
    
    return interactions;
}

std::vector<std::shared_ptr<ComplexFormationInteraction>> MoleculeInteractionLoader::LoadComplexFormationInteractions(const std::string& filePath)
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
                StringDict::stringToId(complexName)
            };
            
            // Create and add interaction (create Molecule objects for proteins)
            interactions.push_back(std::make_shared<ComplexFormationInteraction>(
                Molecule(StringDict::stringToId(firstProtein), ChemicalType::PROTEIN), 
                Molecule(StringDict::stringToId(secondProtein), ChemicalType::PROTEIN), params));
        }
        catch (const std::exception& e) {
            LOG_ERROR("Error parsing complex formation interaction: %s - %s", line.c_str(), e.what());
        }
    }
    
    return interactions;
}

std::vector<std::shared_ptr<TranslationInteraction>> MoleculeInteractionLoader::LoadTranslationInteractions()
{
    std::vector<std::shared_ptr<TranslationInteraction>> interactions;
    
    // Get all mRNA molecules by iterating through StringDict entries
    // We'll create translation interactions for all molecules that have corresponding genes
    for (int i = static_cast<int>(StringDict::ID::eUNKNOWN) + 1; 
         i < static_cast<int>(StringDict::ID::ORGANELLE_END); 
         ++i) {
        
        StringDict::ID id = static_cast<StringDict::ID>(i);
        const std::string& name = StringDict::idToString(id);
        
        // Check if this molecule has GeneData (indicating it can be translated)
        if (GeneWiki::getInstance().hasGeneData(name)) {
            // Create mRNA molecule
            Molecule mRNA(id, ChemicalType::MRNA);
            
            // Get translation rate from MoleculeWiki
            const auto& info = MoleculeWiki::getInfo(mRNA);
            
            // Create translation interaction parameters
            TranslationInteraction::Parameters params{
                info.m_fTranslationRate
            };
            
            // Create and add translation interaction
            interactions.push_back(std::make_shared<TranslationInteraction>(mRNA, params));
        }
    }
    
    return interactions;
}

bool MoleculeInteractionLoader::FileExists(const std::string& filePath)
{
    return std::filesystem::exists(filePath);
}
