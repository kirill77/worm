#include "pch.h"
#include "ProteinWiki.h"
#include "PhosphorylationInteraction.h"
#include "ComplexFormationInteraction.h"
#include "DephosphorylationInteraction.h"
#include <algorithm>
#include <iterator>

// Initialize static member
std::vector<std::shared_ptr<ProteinInteraction>> ProteinWiki::s_proteinInteractions;

void ProteinWiki::Initialize()
{
    // Clear any existing interactions
    s_proteinInteractions.clear();
    
    // === PHOSPHORYLATION INTERACTIONS ===
    
    // PKC-3 (kinase) phosphorylates posterior PARs
    PhosphorylationInteraction::Parameters pkc3ToParParams{
        0.9,    // High removal rate (strong kinase)
        0.07,   // Recovery rate
    };
    
    // PAR-1 (kinase) phosphorylates PAR-3
    PhosphorylationInteraction::Parameters par1ToPar3Params{
        0.7,    // Medium-high removal rate
        0.06,   // Lower recovery rate
    };
    
    // Add phosphorylation interactions
    s_proteinInteractions.push_back(std::make_shared<PhosphorylationInteraction>(
        "PKC-3", "PAR-2", pkc3ToParParams));
    s_proteinInteractions.push_back(std::make_shared<PhosphorylationInteraction>(
        "PKC-3", "PAR-1", pkc3ToParParams));
    s_proteinInteractions.push_back(std::make_shared<PhosphorylationInteraction>(
        "PAR-1", "PAR-3", par1ToPar3Params));
    
    // === DEPHOSPHORYLATION INTERACTIONS ===
    
    // Add dephosphorylation interactions for each protein
    DephosphorylationInteraction::Parameters dephosphoParams{
        0.07,    // Recovery rate
    };
    
    s_proteinInteractions.push_back(std::make_shared<DephosphorylationInteraction>(
        "PAR-2", dephosphoParams));
    s_proteinInteractions.push_back(std::make_shared<DephosphorylationInteraction>(
        "PAR-1", dephosphoParams));
    s_proteinInteractions.push_back(std::make_shared<DephosphorylationInteraction>(
        "PAR-3", dephosphoParams));
    
    // === COMPLEX FORMATION INTERACTIONS ===
    
    // PAR-3 and PAR-6 form a complex
    ComplexFormationInteraction::Parameters par3Par6ComplexParams{
        0.5,             // Binding rate
        0.05,            // Dissociation rate 
        600.0,           // Saturation constant
        "PAR-3-PAR-6"    // Complex name
    };
    
    // PAR-6 and PKC-3 form a complex
    ComplexFormationInteraction::Parameters par6Pkc3ComplexParams{
        0.4,             // Binding rate
        0.04,            // Dissociation rate
        700.0,           // Saturation constant
        "PAR-6-PKC-3"    // Complex name
    };
    
    // Add complex formation interactions
    s_proteinInteractions.push_back(std::make_shared<ComplexFormationInteraction>(
        "PAR-3", "PAR-6", par3Par6ComplexParams));
    s_proteinInteractions.push_back(std::make_shared<ComplexFormationInteraction>(
        "PAR-6", "PKC-3", par6Pkc3ComplexParams));
}

const std::vector<std::shared_ptr<ProteinInteraction>>& ProteinWiki::GetProteinInteractions()
{
    return s_proteinInteractions;
}

std::vector<std::shared_ptr<ProteinInteraction>> ProteinWiki::GetInteractionsInvolvingProtein(const std::string& proteinName)
{
    std::vector<std::shared_ptr<ProteinInteraction>> result;
    
    std::copy_if(s_proteinInteractions.begin(), s_proteinInteractions.end(), 
                 std::back_inserter(result),
                 [&proteinName](const auto& interaction) {
                     return interaction->getProteinA() == proteinName || 
                            interaction->getProteinB() == proteinName;
                 });
                 
    return result;
}

std::vector<std::shared_ptr<ProteinInteraction>> ProteinWiki::GetInteractionsByMechanism(ProteinInteraction::Mechanism mechanism)
{
    std::vector<std::shared_ptr<ProteinInteraction>> result;
    
    std::copy_if(s_proteinInteractions.begin(), s_proteinInteractions.end(), 
                 std::back_inserter(result),
                 [mechanism](const auto& interaction) {
                     return interaction->getMechanism() == mechanism;
                 });
                 
    return result;
} 