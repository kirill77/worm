#pragma once

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include "ProteinInteraction.h"
#include "StringDict.h"
#include "Molecule.h"

// Information about a molecule
struct MolInfo {
    std::string description;      // Description of the molecule
    std::string chemicalFormula;  // Chemical formula (e.g., "C6H12O6")
    double molecularWeight;       // Molecular weight in Daltons
    std::string classification;   // Additional classification info
    double m_fHalfLife;           // How quickly it degrades (in seconds)
    double m_fTranslationRate;    // Rate of protein production
    double m_fChargingRate;       // Rate at which tRNA gets charged with amino acid (for tRNAs only)
    
    MolInfo() : molecularWeight(0.0), m_fHalfLife(0.0), m_fTranslationRate(0.0), m_fChargingRate(0.0) {}
    MolInfo(const std::string& desc, const std::string& formula = "", double weight = 0.0, const std::string& classif = "",
            double halfLife = 0.0, double translationRate = 0.0, double chargingRate = 0.0)
        : description(desc), chemicalFormula(formula), molecularWeight(weight), classification(classif),
          m_fHalfLife(halfLife), m_fTranslationRate(translationRate), m_fChargingRate(chargingRate) {}
};

// A static repository of molecule interaction data
class MoleculeWiki
{
public:
    // Private constructor to prevent instantiation
    MoleculeWiki() = default;
    
private:
    // List of known protein interactions (molecules include proteins)
    static std::vector<std::shared_ptr<ProteinInteraction>> s_proteinInteractions;
    
    // Information database for molecules
    static std::unordered_map<Molecule, MolInfo> m_moleculesInfo;

    // Helper method to load default hardcoded interactions
    static void LoadDefaultInteractions();

public:
    // Initialize all known molecule interactions
    static void Initialize();

    // Get all known protein interactions
    static const std::vector<std::shared_ptr<ProteinInteraction>>& GetProteinInteractions();

    // Get interactions by mechanism
    static std::vector<std::shared_ptr<ProteinInteraction>> GetInteractionsByMechanism(ProteinInteraction::Mechanism mechanism);
    
    // Get the bound molecule name for a molecule on a specific surface
    static std::string GetBoundProteinName(const std::string& proteinName, StringDict::ID surface);
    
    // Get information about a specific molecule
    static const MolInfo& getInfo(const Molecule& molecule);
    
    // Initialize tRNA molecule information with charging rates
    static void initializeTRNAInfo();
};
