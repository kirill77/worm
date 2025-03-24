#include "pch.h"
#include "ProteinWiki.h"
#include "PhosphorylationInteraction.h"
#include "ComplexFormationInteraction.h"
#include "DephosphorylationInteraction.h"
#include "ProteinBindingSurface.h"
#include <algorithm>
#include <iterator>

// Initialize static members
std::vector<std::shared_ptr<ProteinInteraction>> ProteinWiki::s_proteinInteractions;

void ProteinWiki::Initialize()
{
    // Clear any existing interactions
    s_proteinInteractions.clear();
    
    // === PHOSPHORYLATION INTERACTIONS ===
    
    // PKC-3 (kinase) phosphorylates posterior PARs, but only when in complex with PAR-6
    PhosphorylationInteraction::Parameters pkc3ComplexToParParams{
        0.9,    // High removal rate (strong kinase)
        0.07,   // Saturation constant for Hill-type kinetics
    };
    
    // PAR-1 (kinase) phosphorylates PAR-3
    PhosphorylationInteraction::Parameters par1ToPar3Params{
        0.7,    // Medium-high removal rate
        0.06,   // Saturation constant for Hill-type kinetics
    };
    
    // Add phosphorylation interactions - using the PAR-6-PKC-3 complex as the active kinase
    s_proteinInteractions.push_back(std::make_shared<PhosphorylationInteraction>(
        "PAR-6-PKC-3", "PAR-2", pkc3ComplexToParParams));
    s_proteinInteractions.push_back(std::make_shared<PhosphorylationInteraction>(
        "PAR-6-PKC-3", "PAR-1", pkc3ComplexToParParams));
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
          
    // Define binding site name for membrane
    std::string bindingSiteName = GetBindingSiteName(BindingSurface::MEMBRANE);
    
    // Parameters for PAR protein binding to membrane (converted to ComplexFormationInteraction parameters)
    ComplexFormationInteraction::Parameters par1MembraneParams{
        0.6,    // High binding rate for posterior proteins
        0.04,   // Moderate dissociation rate
        800.0,  // Saturation constant
        GetBoundProteinName("PAR-1", BindingSurface::MEMBRANE)  // Complex name
    };
    
    ComplexFormationInteraction::Parameters par2MembraneParams{
        0.5,    // Medium binding rate
        0.03,   // Low dissociation rate
        900.0,  // Saturation constant
        GetBoundProteinName("PAR-2", BindingSurface::MEMBRANE)  // Complex name
    };
    
    ComplexFormationInteraction::Parameters par3MembraneParams{
        0.4,    // Lower binding rate for anterior proteins
        0.1,    // Higher dissociation rate
        1000.0, // Saturation constant
        GetBoundProteinName("PAR-3", BindingSurface::MEMBRANE)  // Complex name
    };
    
    // Add membrane binding interactions (now using ComplexFormationInteraction)
    s_proteinInteractions.push_back(std::make_shared<ComplexFormationInteraction>(
        "PAR-1", bindingSiteName, par1MembraneParams));
    s_proteinInteractions.push_back(std::make_shared<ComplexFormationInteraction>(
        "PAR-2", bindingSiteName, par2MembraneParams));
    s_proteinInteractions.push_back(std::make_shared<ComplexFormationInteraction>(
        "PAR-3", bindingSiteName, par3MembraneParams));
}

const std::vector<std::shared_ptr<ProteinInteraction>>& ProteinWiki::GetProteinInteractions()
{
    return s_proteinInteractions;
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

std::string ProteinWiki::GetPhosphorylatedName(const std::string& proteinName)
{
    return proteinName + "-P";
}

std::string ProteinWiki::BindingSurfaceToString(BindingSurface surface)
{
    switch (surface)
    {
        case BindingSurface::MEMBRANE:
            return "MEMBRANE";
        case BindingSurface::CORTEX:
            return "CORTEX";
        case BindingSurface::CENTROSOME:
            return "CENTROSOME";
        default:
            return "UNKNOWN";
    }
}

std::string ProteinWiki::GetBindingSiteName(BindingSurface surface)
{
    return "BINDING-SITE-" + BindingSurfaceToString(surface);
}

std::string ProteinWiki::GetBoundProteinName(const std::string& proteinName, BindingSurface surface)
{
    return proteinName + "-" + BindingSurfaceToString(surface);
}
