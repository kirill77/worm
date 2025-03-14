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
    , m_lastPhosphorylatedAmount(0.0)
{
}

bool PhosphorylationInteraction::apply(GridCell& cell, double dt, double& atpConsumed) const
{
    // Reset tracking variable
    m_lastPhosphorylatedAmount = 0.0;
    
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
        
        // Remove phosphorylated proteins from active population
        targetIt->second.m_fNumber -= phosphorylatedAmount;
        
        // Track how much was phosphorylated for recovery effects
        m_lastPhosphorylatedAmount = phosphorylatedAmount;
        
        return true;
    }
    
    return false;
}

bool PhosphorylationInteraction::applyNeighborEffects(
    GridCell& cell, 
    std::vector<std::reference_wrapper<GridCell>>& neighborCells, 
    double dt) const
{
    // If no phosphorylation occurred, no recovery needed
    if (m_lastPhosphorylatedAmount <= 0 || neighborCells.empty()) {
        return false;
    }
    
    // Calculate recovery
    double recoveredAmount = m_lastPhosphorylatedAmount * m_recoveryRate * dt;
    
    if (recoveredAmount <= 0) {
        return false;
    }
    
    // Count non-cortical neighbors (assuming we have a function to check this)
    int nonCorticalCount = 0;
    for (auto& neighborRef : neighborCells) {
        // For this example, let's assume we're distributing to all neighbors
        // In your real implementation, you can filter based on cortical status
        nonCorticalCount++;
    }
    
    if (nonCorticalCount == 0) {
        return false;
    }
    
    // Distribute recovery to neighbors
    double amountPerNeighbor = recoveredAmount / nonCorticalCount;
    
    for (auto& neighborRef : neighborCells) {
        GridCell& neighbor = neighborRef.get();
        
        // Add recovered proteins to neighbor
        auto& neighborPop = neighbor.getOrCreateProtein(m_proteinB);
        neighborPop.m_fNumber += amountPerNeighbor;
    }
    
    return true;
} 