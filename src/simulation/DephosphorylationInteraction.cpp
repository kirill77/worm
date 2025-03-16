#include "pch.h"
#include "DephosphorylationInteraction.h"
#include <algorithm>

DephosphorylationInteraction::DephosphorylationInteraction(
    const std::string& target,
    const Parameters& params)
    : ProteinInteraction("", target,  // No kinase needed for dephosphorylation
                        Mechanism::DEPHOSPHORYLATION,
                        0.1)  // Lower ATP cost for dephosphorylation
    , m_recoveryRate(params.recoveryRate)
{
}

bool DephosphorylationInteraction::apply(GridCell& cell, double dt, double& atpConsumed) const
{
    // Get phosphorylated protein population (e.g. "PAR-2-P")
    std::string phosphorylatedName = m_proteinB + "-P";
    auto phosphorylatedIt = cell.m_proteins.find(phosphorylatedName);
    
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
    
    // Check ATP availability
    if (cell.m_fAtp < requiredATP) {
        // Scale down recovery based on available ATP
        recoveredAmount = recoveredAmount * (cell.m_fAtp / requiredATP);
        requiredATP = cell.m_fAtp;
    }
    
    // Remove from phosphorylated population
    phosphorylatedIt->second.m_fNumber -= recoveredAmount;
    
    // Add back to original unphosphorylated population
    auto& unphosphorylatedPop = cell.getOrCreateProtein(m_proteinB);
    unphosphorylatedPop.m_fNumber += recoveredAmount;
    
    // Update ATP consumption
    atpConsumed += requiredATP;
    cell.m_fAtp -= requiredATP;
    
    return true;
} 