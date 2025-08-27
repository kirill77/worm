#include "pch.h"
#include "PhosphorylationInteraction.h"
#include "MoleculeWiki.h"
#include "ResourceDistributor.h"
#include <algorithm>
#include <cmath>

PhosphorylationInteraction::PhosphorylationInteraction(
    const std::string& kinase, 
    const std::string& target, 
    const Parameters& params)
    : ProteinInteraction(Mechanism::PHOSPHORYLATION, 0.5)  // Standard ATP cost for phosphorylation
    , m_kinaseName(kinase)
    , m_targetName(target)
    , m_phosphorylatedName(MoleculeWiki::GetPhosphorylatedName(target))
    , m_removalRate(params.removalRate)
    , m_saturationConstant(params.saturationConstant)
{
}

bool PhosphorylationInteraction::apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const
{
    double kinaseAmount = resDistributor.getAvailableResource(m_kinaseName);
    double targetAmount = resDistributor.getAvailableResource(m_targetName);
    
    // Calculate phosphorylation using Hill-like kinetics
    double removalRate = m_removalRate * kinaseAmount / (m_saturationConstant + kinaseAmount);
    
    // Calculate amount to phosphorylate in this time step
    double phosphorylatedAmount = removalRate * targetAmount * dt;
    
    // Phosphorylation requires ATP
    double requiredATP = phosphorylatedAmount * m_atpCost;
    
    // If we're in a dry run, just report resource requirements and return
    if (resDistributor.isDryRun()) {
        if (phosphorylatedAmount > 0) {
            // Register our requirements with the resource distributor
            resDistributor.notifyResourceWanted("ATP", requiredATP);
            resDistributor.notifyResourceWanted(m_targetName, phosphorylatedAmount);
            return true; // We're reporting resource needs
        }
        return false;
    }

    // Apply the effect if any phosphorylation occurs
    if (phosphorylatedAmount > 0) {
        // Update ATP consumption
        auto& atpMolecule = cell.getOrCreateMolecule("ATP");
        atpMolecule.m_fNumber -= requiredATP;
        assert(atpMolecule.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert ATP doesn't go below minimum
        
        // Remove proteins from unphosphorylated population
        auto targetIt = cell.m_molecules.find(m_targetName);
        targetIt->second.m_fNumber -= phosphorylatedAmount;
        assert(targetIt->second.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert protein level doesn't go below minimum
        
        // Add to phosphorylated population
        auto& phosphorylatedPop = cell.getOrCreateMolecule(m_phosphorylatedName);
        phosphorylatedPop.m_fNumber += phosphorylatedAmount;
        
        return true;
    }
    
    return false;
}
