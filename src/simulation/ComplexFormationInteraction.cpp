#include "pch.h"
#include "ComplexFormationInteraction.h"
#include "ResourceDistributor.h"
#include <algorithm>
#include <cmath>

ComplexFormationInteraction::ComplexFormationInteraction(
    const std::string& firstProtein, 
    const std::string& secondProtein, 
    const Parameters& params)
    : ProteinInteraction(Mechanism::BINDING, 0.2)  // Lower ATP cost for binding
    , m_firstProteinName(firstProtein)
    , m_secondProteinName(secondProtein)
    , m_bindingRate(params.bindingRate)
    , m_dissociationRate(params.dissociationRate)
    , m_saturationConstant(params.saturationConstant)
    , m_complexName(params.complexName)
{
}

bool ComplexFormationInteraction::apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const
{
    // Check for both proteins
    auto firstProteinIt = cell.m_proteins.find(m_firstProteinName);
    auto secondProteinIt = cell.m_proteins.find(m_secondProteinName);
    
    if (firstProteinIt == cell.m_proteins.end() || firstProteinIt->second.m_fNumber <= 0 ||
        secondProteinIt == cell.m_proteins.end() || secondProteinIt->second.m_fNumber <= 0) {
        return false; // One or both proteins missing
    }
    
    double firstProteinAmount = resDistributor.getAvailableResource(m_firstProteinName);
    double secondProteinAmount = resDistributor.getAvailableResource(m_secondProteinName);
    
    // Calculate binding using mass action kinetics
    double bindingPotential = m_bindingRate * firstProteinAmount * secondProteinAmount / 
                             (m_saturationConstant + firstProteinAmount + secondProteinAmount);
    
    // The amount that can actually bind is limited by the lesser of the two proteins
    double boundAmount = std::min(bindingPotential * dt, std::min(firstProteinAmount, secondProteinAmount));
    
    // Binding requires ATP
    double requiredATP = boundAmount * m_atpCost;
    
    // Also check for dissociation of existing complexes
    auto complexIt = cell.m_proteins.find(m_complexName);
    double complexAmount = (complexIt != cell.m_proteins.end()) ? complexIt->second.m_fNumber : 0.0;
    
    // Calculate dissociation (simpler first-order kinetics)
    double dissociatedAmount = complexAmount * m_dissociationRate * dt;
    
    // If we're in a dry run, just report resource requirements and return
    if (resDistributor.isDryRun()) {
        if (boundAmount > 0) {
            // Register our requirements with the resource distributor
            resDistributor.notifyResourceWanted("ATP", requiredATP);
            resDistributor.notifyResourceWanted(m_firstProteinName, boundAmount);
            resDistributor.notifyResourceWanted(m_secondProteinName, boundAmount);
            return true; // We're reporting resource needs
        }
        // Dissociation doesn't consume resources, but still return true if it occurs
        return dissociatedAmount > 0;
    }

    // Apply binding if any occurs
    if (boundAmount > 0) {
        // Update ATP consumption
        cell.m_fAtp -= requiredATP;
        
        // Remove proteins from free populations
        firstProteinIt->second.m_fNumber -= boundAmount;
        secondProteinIt->second.m_fNumber -= boundAmount;
        
        // Add to complex population
        auto& complexPop = cell.getOrCreateProtein(m_complexName);
        complexPop.m_fNumber += boundAmount;
    }
    
    // Apply dissociation if any occurs
    if (dissociatedAmount > 0 && complexIt != cell.m_proteins.end()) {
        // Remove from complex population
        complexIt->second.m_fNumber -= dissociatedAmount;
        
        // Return to free protein populations
        firstProteinIt->second.m_fNumber += dissociatedAmount;
        secondProteinIt->second.m_fNumber += dissociatedAmount;
    }
    
    return (boundAmount > 0 || dissociatedAmount > 0);
}
