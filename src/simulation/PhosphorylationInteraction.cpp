#include "pch.h"
#include "PhosphorylationInteraction.h"
#include <algorithm>
#include <cmath>

PhosphorylationInteraction::PhosphorylationInteraction(
    const std::string& kinase, 
    const std::string& target, 
    const Parameters& params)
    : ProteinInteraction(kinase, target, 
                       Mechanism::PHOSPHORYLATION, 
                       0.5)  // Standard ATP cost for phosphorylation
    , m_removalRate(params.removalRate)
    , m_recoveryRate(params.recoveryRate)
    , m_saturationConstant(params.saturationConstant)
{
}

bool PhosphorylationInteraction::apply(GridCell& cell, double dt, double& atpConsumed) const
{
    // Get kinase amount
    auto kinaseIt = cell.m_proteins.find(m_proteinA);
    if (kinaseIt == cell.m_proteins.end() || kinaseIt->second.m_fNumber <= 0) {
        return false; // No kinase present
    }
    
    // Get target amount
    auto targetIt = cell.m_proteins.find(m_proteinB);
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
    
    // Check ATP availability
    if (cell.m_fAtp < requiredATP) {
        // Scale down phosphorylation based on available ATP
        phosphorylatedAmount = phosphorylatedAmount * (cell.m_fAtp / requiredATP);
        requiredATP = cell.m_fAtp;
    }
    
    // Apply the effect if any phosphorylation occurs
    if (phosphorylatedAmount > 0) {
        // Update ATP consumption
        atpConsumed += requiredATP;
        cell.m_fAtp -= requiredATP;
        
        // Remove proteins from unphosphorylated population
        targetIt->second.m_fNumber -= phosphorylatedAmount;
        
        // Add to phosphorylated population
        std::string phosphorylatedName = m_proteinB + "-P";  // e.g., "PAR-2" becomes "PAR-2-P"
        auto& phosphorylatedPop = cell.getOrCreateProtein(phosphorylatedName);
        phosphorylatedPop.m_fNumber += phosphorylatedAmount;
        
        return true;
    }
    
    return false;
}

bool PhosphorylationInteraction::applyNeighborEffects(
    GridCell& cell, 
    std::vector<std::reference_wrapper<GridCell>>& neighborCells, 
    double dt) const
{
    // Get phosphorylated protein population
    std::string phosphorylatedName = m_proteinB + "-P";
    auto phosphorylatedIt = cell.m_proteins.find(phosphorylatedName);
    
    if (phosphorylatedIt == cell.m_proteins.end() || 
        phosphorylatedIt->second.m_fNumber <= 0 || 
        neighborCells.empty()) {
        return false;
    }
    
    // Calculate recovery
    double phosphorylatedAmount = phosphorylatedIt->second.m_fNumber;
    double recoveredAmount = phosphorylatedAmount * m_recoveryRate * dt;
    
    if (recoveredAmount <= 0) {
        return false;
    }
    
    // Remove from phosphorylated population
    phosphorylatedIt->second.m_fNumber -= recoveredAmount;
    
    // Distribute recovered (unphosphorylated) proteins to neighbors
    double amountPerNeighbor = recoveredAmount / neighborCells.size();
    
    for (auto& neighborRef : neighborCells) {
        GridCell& neighbor = neighborRef.get();
        // Add back to original unphosphorylated protein population
        auto& neighborPop = neighbor.getOrCreateProtein(m_proteinB);
        neighborPop.m_fNumber += amountPerNeighbor;
    }
    
    return true;
} 