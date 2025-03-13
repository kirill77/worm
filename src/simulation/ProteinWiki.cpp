#include "pch.h"
#include "ProteinWiki.h"
#include <algorithm>
#include <iterator>

// Initialize static member
std::vector<ProteinAntagonism> ProteinWiki::s_proteinAntagonisms;

void ProteinWiki::Initialize()
{
    // Clear any existing antagonisms
    s_proteinAntagonisms.clear();
    
    // PKC-3 (kinase) phosphorylates posterior PARs
    ProteinAntagonism::Parameters pkc3ToParParams{
        0.9,                                          // High removal rate (strong kinase)
        0.07,                                         // Recovery rate
        550.0,                                        // Low saturation for stronger effect
        ProteinAntagonism::Mechanism::PHOSPHORYLATION, // Phosphorylation mechanism
        0.5                                           // ATP cost
    };
    
    // PAR-1 (kinase) phosphorylates PAR-3
    ProteinAntagonism::Parameters par1ToPar3Params{
        0.7,                                          // Medium-high removal rate
        0.06,                                         // Lower recovery rate
        650.0,                                        // Medium saturation constant
        ProteinAntagonism::Mechanism::PHOSPHORYLATION, // Phosphorylation mechanism
        0.4                                           // ATP cost
    };
    
    // PAR-1 weakly affects PAR-6 (indirect)
    ProteinAntagonism::Parameters par1ToPar6Params{
        0.3,                                          // Lower removal rate
        0.1,                                          // Higher recovery rate
        900.0,                                        // High saturation constant
        ProteinAntagonism::Mechanism::RECRUITMENT,     // Indirect mechanism
        0.0                                           // No ATP cost
    };
    
    // PAR-2 affects anterior proteins through cortical exclusion
    ProteinAntagonism::Parameters par2ToPar3Params{
        0.5,                                          // Medium removal rate 
        0.09,                                         // Medium recovery rate
        750.0,                                        // Medium saturation constant
        ProteinAntagonism::Mechanism::CORTICAL_EXCLUSION, // Cortical competition
        0.0                                           // No ATP cost
    };
    
    ProteinAntagonism::Parameters par2ToPar6Params{
        0.35,                                         // Medium-low removal rate
        0.09,                                         // Recovery rate
        800.0,                                        // Medium-high saturation constant
        ProteinAntagonism::Mechanism::CORTICAL_EXCLUSION, // Cortical competition
        0.0                                           // No ATP cost
    };
    
    ProteinAntagonism::Parameters par2ToPkc3Params{
        0.3,                                          // Lower removal rate
        0.1,                                          // Higher recovery rate
        850.0,                                        // Higher saturation constant
        ProteinAntagonism::Mechanism::CORTICAL_EXCLUSION, // Cortical competition
        0.0                                           // No ATP cost
    };
    
    // Add antagonistic relationships
    s_proteinAntagonisms.emplace_back("PKC-3", "PAR-2", pkc3ToParParams);
    s_proteinAntagonisms.emplace_back("PKC-3", "PAR-1", pkc3ToParParams);
    
    s_proteinAntagonisms.emplace_back("PAR-1", "PAR-3", par1ToPar3Params);
    s_proteinAntagonisms.emplace_back("PAR-1", "PAR-6", par1ToPar6Params);
    
    s_proteinAntagonisms.emplace_back("PAR-2", "PAR-3", par2ToPar3Params);
    s_proteinAntagonisms.emplace_back("PAR-2", "PAR-6", par2ToPar6Params);
    s_proteinAntagonisms.emplace_back("PAR-2", "PKC-3", par2ToPkc3Params);
}

const std::vector<ProteinAntagonism>& ProteinWiki::GetProteinAntagonisms()
{
    return s_proteinAntagonisms;
}

std::vector<ProteinAntagonism> ProteinWiki::GetAntagonismsFromProtein(const std::string& antagonistName)
{
    std::vector<ProteinAntagonism> result;
    
    std::copy_if(s_proteinAntagonisms.begin(), s_proteinAntagonisms.end(), 
                 std::back_inserter(result),
                 [&antagonistName](const ProteinAntagonism& antagonism) {
                     return antagonism.getAntagonist() == antagonistName;
                 });
                 
    return result;
}

std::vector<ProteinAntagonism> ProteinWiki::GetAntagonismsToProtein(const std::string& targetName)
{
    std::vector<ProteinAntagonism> result;
    
    std::copy_if(s_proteinAntagonisms.begin(), s_proteinAntagonisms.end(), 
                 std::back_inserter(result),
                 [&targetName](const ProteinAntagonism& antagonism) {
                     return antagonism.getTarget() == targetName;
                 });
                 
    return result;
} 