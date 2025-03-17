#include "pch.h"
#include "SurfaceBindingInteraction.h"
#include <algorithm>
#include <cmath>

SurfaceBindingInteraction::SurfaceBindingInteraction(
    const std::string& proteinName, 
    const std::string& bindingSiteName,
    const std::string& boundComplexName,
    const Parameters& params)
    : ProteinInteraction(Mechanism::BINDING, 0.1)  // Low ATP cost for binding to surface
    , m_proteinName(proteinName)
    , m_bindingSiteName(bindingSiteName)
    , m_boundComplexName(boundComplexName)
    , m_bindingRate(params.bindingRate)
    , m_dissociationRate(params.dissociationRate)
    , m_saturationConstant(params.saturationConstant)
{
}

bool SurfaceBindingInteraction::apply(GridCell& cell, double dt, double& atpConsumed) const
{
    // Check for protein and binding sites
    auto proteinIt = cell.m_proteins.find(m_proteinName);
    auto bindingSiteIt = cell.m_proteins.find(m_bindingSiteName);
    
    if (proteinIt == cell.m_proteins.end() || proteinIt->second.m_fNumber <= 0 ||
        bindingSiteIt == cell.m_proteins.end() || bindingSiteIt->second.m_fNumber <= 0) {
        return false; // Missing protein or binding sites
    }
    
    double proteinAmount = proteinIt->second.m_fNumber;
    double bindingSiteAmount = bindingSiteIt->second.m_fNumber;
    
    // Get current amount of bound complex, if any
    auto boundIt = cell.m_proteins.find(m_boundComplexName);
    double boundAmount = (boundIt != cell.m_proteins.end()) ? boundIt->second.m_fNumber : 0.0;
    
    // Calculate binding using mass action kinetics with saturation
    double bindingPotential = m_bindingRate * proteinAmount * bindingSiteAmount / 
                             (m_saturationConstant + proteinAmount);
    
    // Calculate amount to bind in this time step (limited by available protein and binding sites)
    double newBoundAmount = std::min(bindingPotential * dt, std::min(proteinAmount, bindingSiteAmount));
    
    // Calculate dissociation of existing bound complexes
    double dissociatedAmount = boundAmount * m_dissociationRate * dt;
    
    // Binding requires ATP
    double requiredATP = newBoundAmount * m_atpCost;
    
    // Check ATP availability
    if (cell.m_fAtp < requiredATP) {
        // Scale down binding based on available ATP
        newBoundAmount = newBoundAmount * (cell.m_fAtp / requiredATP);
        requiredATP = cell.m_fAtp;
    }
    
    bool changesApplied = false;
    
    // Apply binding if any occurs
    if (newBoundAmount > 0) {
        // Update ATP consumption
        atpConsumed += requiredATP;
        cell.m_fAtp -= requiredATP;
        
        // Remove proteins from free populations
        proteinIt->second.m_fNumber -= newBoundAmount;
        bindingSiteIt->second.m_fNumber -= newBoundAmount;  // Occupy binding sites
        
        // Add to bound complex population
        auto& boundPop = cell.getOrCreateProtein(m_boundComplexName);
        // bind this protein to the surface containing the binding site
        boundPop.bindTo(bindingSiteIt->second.getBindingSurface());
        boundPop.m_fNumber += newBoundAmount;
        
        changesApplied = true;
    }
    
    // Apply dissociation if any occurs
    if (dissociatedAmount > 0 && boundIt != cell.m_proteins.end()) {
        // Remove from bound complex population
        boundIt->second.m_fNumber -= dissociatedAmount;
        
        // Return to free populations
        proteinIt->second.m_fNumber += dissociatedAmount;
        bindingSiteIt->second.m_fNumber += dissociatedAmount;  // Free up binding sites
        
        changesApplied = true;
    }
    
    return changesApplied;
} 