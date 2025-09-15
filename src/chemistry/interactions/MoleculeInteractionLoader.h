#pragma once

#include "PhosphorylationInteraction.h"
#include "DephosphorylationInteraction.h"
#include "ComplexFormationInteraction.h"
#include "TranslationInteraction.h"
#include "chemistry/molecules/MoleculeWiki.h"
#include <string>
#include <vector>
#include <memory>
#include <filesystem>

/**
 * @class MoleculeInteractionLoader
 * @brief Loads molecule interactions from CSV files
 * 
 * This class is responsible for loading various types of molecule interactions
 * from their respective CSV configuration files.
 */
class MoleculeInteractionLoader
{
public:
    /**
     * @brief Load all molecule interactions from CSV files
     * @param basePath Base directory containing the CSV files
     * @return Vector of loaded molecule interactions
     */
    static std::vector<std::shared_ptr<MoleculeInteraction>> LoadAllInteractions(
        const std::string& basePath = "data/interactions/");

private:
    /**
     * @brief Load phosphorylation interactions from CSV
     * @param filePath Path to the CSV file
     * @return Vector of loaded phosphorylation interactions
     */
    static std::vector<std::shared_ptr<PhosphorylationInteraction>> LoadPhosphorylationInteractions(
        const std::string& filePath);
    
    /**
     * @brief Load dephosphorylation interactions from CSV
     * @param filePath Path to the CSV file
     * @return Vector of loaded dephosphorylation interactions
     */
    static std::vector<std::shared_ptr<DephosphorylationInteraction>> LoadDephosphorylationInteractions(
        const std::string& filePath);
    
    /**
     * @brief Load complex formation interactions from CSV
     * @param filePath Path to the CSV file
     * @return Vector of loaded complex formation interactions
     */
    static std::vector<std::shared_ptr<ComplexFormationInteraction>> LoadComplexFormationInteractions(
        const std::string& filePath);
    
    /**
     * @brief Load translation interactions for all mRNA molecules
     * @return Vector of loaded translation interactions
     */
    static std::vector<std::shared_ptr<TranslationInteraction>> LoadTranslationInteractions();
    
    /**
     * @brief Check if the file exists
     * @param filePath Path to the file to check
     * @return True if the file exists, false otherwise
     */
    static bool FileExists(const std::string& filePath);
};
