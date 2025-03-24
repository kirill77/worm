#include "pch.h"
#include "PhosphorylationInteraction.h"
#include "ProteinWiki.h"
#include "ResourceAllocation.h"
#include <algorithm>
#include <cmath>

PhosphorylationInteraction::PhosphorylationInteraction(
    const std::string& kinase, 
    const std::string& target, 
    const Parameters& params)
    : ProteinInteraction(Mechanism::PHOSPHORYLATION, 0.5)  // Standard ATP cost for phosphorylation
    , m_kinaseName(kinase)
    , m_targetName(target)
    , m_phosphorylatedName(ProteinWiki::GetPhosphorylatedName(target))
    , m_removalRate(params.removalRate)
    , m_saturationConstant(params.saturationConstant)
{
}

bool PhosphorylationInteraction::apply(GridCell& cell, double dt, ResourceAllocation& resDistributor) const
{
    // Get kinase amount
    auto kinaseIt = cell.m_proteins.find(m_kinaseName);
    if (kinaseIt == cell.m_proteins.end() || kinaseIt->second.m_fNumber <= 0) {
        return false; // No kinase present
    }
    
    // Get target amount
    auto targetIt = cell.m_proteins.find(m_targetName);
    if (targetIt == cell.m_proteins.end() || targetIt->second.m_fNumber <= 0) {
        return false; // No target present
    }
    
    double kinaseAmount = kinaseIt->second.m_fNumber;
    double targetAmount = targetIt->second.m_fNumber;
    
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
            resDistributor.notifyResourceConsumed("ATP", requiredATP);
            resDistributor.notifyResourceConsumed(m_targetName, phosphorylatedAmount);
            return true; // We're reporting resource needs
        }
        return false;
    }
    
    // This is the real run, get the scaling factor for this interaction
    double scalingFactor = resDistributor.notifyNewInteractionStarting(*this);
    
    // Scale our resource usage by the scaling factor
    phosphorylatedAmount *= scalingFactor;
    requiredATP *= scalingFactor;
    
    // Apply the effect if any phosphorylation occurs
    if (phosphorylatedAmount > 0) {
        // Update ATP consumption
        cell.m_fAtp -= requiredATP;
        
        // Remove proteins from unphosphorylated population
        targetIt->second.m_fNumber -= phosphorylatedAmount;
        
        // Add to phosphorylated population
        auto& phosphorylatedPop = cell.getOrCreateProtein(m_phosphorylatedName);
        phosphorylatedPop.m_fNumber += phosphorylatedAmount;
        
        return true;
    }
    
    return false;
}
