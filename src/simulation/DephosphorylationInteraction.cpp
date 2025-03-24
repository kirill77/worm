#include "pch.h"
#include "DephosphorylationInteraction.h"
#include "ProteinWiki.h"
#include "ResourceAllocation.h"
#include <algorithm>

DephosphorylationInteraction::DephosphorylationInteraction(
    const std::string& target,
    const Parameters& params)
    : ProteinInteraction(Mechanism::DEPHOSPHORYLATION, 0.1)  // Lower ATP cost for dephosphorylation
    , m_targetName(target)
    , m_phosphorylatedName(ProteinWiki::GetPhosphorylatedName(target))
    , m_recoveryRate(params.recoveryRate)
{
}

bool DephosphorylationInteraction::apply(GridCell& cell, double dt, ResourceAllocation& resDistributor) const
{
    // Get phosphorylated protein population
    auto phosphorylatedIt = cell.m_proteins.find(m_phosphorylatedName);
    
    if (phosphorylatedIt == cell.m_proteins.end() || 
        phosphorylatedIt->second.m_fNumber <= 0) {
        return false;
    }
    
    // Calculate recovery
    double phosphorylatedAmount = phosphorylatedIt->second.m_fNumber;
    double recoveredAmount = phosphorylatedAmount * m_recoveryRate * dt;
    
    if (recoveredAmount <= 0) {
        return false;
    }
    
    // Dephosphorylation requires a small amount of ATP
    double requiredATP = recoveredAmount * m_atpCost;
    
    // If we're in a dry run, just report resource requirements and return
    if (resDistributor.isDryRun()) {
        // Register our requirements with the resource distributor
        resDistributor.notifyResourceConsumed("ATP", requiredATP);
        resDistributor.notifyResourceConsumed(m_phosphorylatedName, recoveredAmount);
        return true; // We're reporting resource needs
    }
    
    // This is the real run, get the scaling factor for this interaction
    double scalingFactor = resDistributor.notifyNewInteractionStarting(*this);
    
    // Scale our resource usage by the scaling factor
    recoveredAmount *= scalingFactor;
    requiredATP *= scalingFactor;
    
    // Remove from phosphorylated population
    phosphorylatedIt->second.m_fNumber -= recoveredAmount;
    
    // Add back to original unphosphorylated population
    auto& unphosphorylatedPop = cell.getOrCreateProtein(m_targetName);
    unphosphorylatedPop.m_fNumber += recoveredAmount;
    
    // Update ATP consumption
    cell.m_fAtp -= requiredATP;
    
    return true;
} 