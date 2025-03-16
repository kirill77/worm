#include "pch.h"
#include "ComplexFormationInteraction.h"
#include <algorithm>
#include <cmath>

ComplexFormationInteraction::ComplexFormationInteraction(
    const std::string& proteinA, 
    const std::string& proteinB, 
    const Parameters& params)
    : ProteinInteraction(proteinA, proteinB, 
                       Mechanism::BINDING, 
                       0.2)  // Lower ATP cost for binding
    , m_bindingRate(params.bindingRate)
    , m_dissociationRate(params.dissociationRate)
    , m_saturationConstant(params.saturationConstant)
    , m_complexName(params.complexName)
{
}

bool ComplexFormationInteraction::apply(GridCell& cell, double dt, double& atpConsumed) const
{
    // Check for both proteins
    auto proteinAIt = cell.m_proteins.find(m_proteinA);
    auto proteinBIt = cell.m_proteins.find(m_proteinB);
    
    if (proteinAIt == cell.m_proteins.end() || proteinAIt->second.m_fNumber <= 0 ||
        proteinBIt == cell.m_proteins.end() || proteinBIt->second.m_fNumber <= 0) {
        return false; // One or both proteins missing
    }
    
    double proteinAAmount = proteinAIt->second.m_fNumber;
    double proteinBAmount = proteinBIt->second.m_fNumber;
    
    // Calculate binding using mass action kinetics
    double bindingPotential = m_bindingRate * proteinAAmount * proteinBAmount / 
                             (m_saturationConstant + proteinAAmount + proteinBAmount);
    
    // The amount that can actually bind is limited by the lesser of the two proteins
    double boundAmount = std::min(bindingPotential * dt, std::min(proteinAAmount, proteinBAmount));
    
    // Binding requires ATP
    double requiredATP = boundAmount * m_atpCost;
    
    // Check ATP availability
    if (cell.m_fAtp < requiredATP) {
        // Scale down binding based on available ATP
        boundAmount = boundAmount * (cell.m_fAtp / requiredATP);
        requiredATP = cell.m_fAtp;
    }
    
    // Also check for dissociation of existing complexes
    auto complexIt = cell.m_proteins.find(m_complexName);
    double complexAmount = (complexIt != cell.m_proteins.end()) ? complexIt->second.m_fNumber : 0.0;
    
    // Calculate dissociation (simpler first-order kinetics)
    double dissociatedAmount = complexAmount * m_dissociationRate * dt;
    
    // Apply binding if any occurs
    if (boundAmount > 0) {
        // Update ATP consumption
        atpConsumed += requiredATP;
        cell.m_fAtp -= requiredATP;
        
        // Remove proteins from free populations
        proteinAIt->second.m_fNumber -= boundAmount;
        proteinBIt->second.m_fNumber -= boundAmount;
        
        // Add to complex population
        auto& complexPop = cell.getOrCreateProtein(m_complexName);
        complexPop.m_fNumber += boundAmount;
    }
    
    // Apply dissociation if any occurs
    if (dissociatedAmount > 0 && complexIt != cell.m_proteins.end()) {
        // Remove from complex population
        complexIt->second.m_fNumber -= dissociatedAmount;
        
        // Return to free protein populations
        proteinAIt->second.m_fNumber += dissociatedAmount;
        proteinBIt->second.m_fNumber += dissociatedAmount;
    }
    
    return (boundAmount > 0 || dissociatedAmount > 0);
}
