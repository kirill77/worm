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
    auto complexIt = cell.m_molecules.find(Molecule(m_complexName, ChemicalType::PROTEIN));
    double complexAmount = (complexIt != cell.m_molecules.end()) ? complexIt->second.m_fNumber : 0.0;
    
    // Calculate dissociation (simpler first-order kinetics)
    double dissociatedAmount = complexAmount * m_dissociationRate * dt;
    
    // If we're in a dry run, just report resource requirements and return
    if (resDistributor.isDryRun()) {
        if (boundAmount > 0) {
            // Register our requirements with the resource distributor
            resDistributor.notifyResourceWanted(StringDict::idToString(StringDict::ID::ATP), requiredATP);
            resDistributor.notifyResourceWanted(m_firstProteinName, boundAmount);
            resDistributor.notifyResourceWanted(m_secondProteinName, boundAmount);
            return true; // We're reporting resource needs
        }
        // Dissociation doesn't consume resources, but still return true if it occurs
        return dissociatedAmount > 0;
    }

    auto firstProteinIt = cell.m_molecules.find(Molecule(m_firstProteinName, ChemicalType::PROTEIN));
    auto secondProteinIt = cell.m_molecules.find(Molecule(m_secondProteinName, ChemicalType::PROTEIN));

    // Apply binding if any occurs
    if (boundAmount > 0) {
        // Update ATP consumption
        auto& atpPop = cell.getOrCreateMolPop(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE));
        atpPop.m_fNumber -= requiredATP;
        assert(atpPop.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert ATP doesn't go below minimum
        
        // Remove proteins from free populations
        firstProteinIt->second.m_fNumber -= boundAmount;
        assert(firstProteinIt->second.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert protein level doesn't go below minimum
        
        secondProteinIt->second.m_fNumber -= boundAmount;
        assert(secondProteinIt->second.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert protein level doesn't go below minimum
        
        // Add to complex population
        auto& complexPop = cell.getOrCreateMolPop(Molecule(m_complexName, ChemicalType::PROTEIN));
        complexPop.m_fNumber += boundAmount;

        // update binding surface
        assert(!firstProteinIt->second.isBound());
        if (secondProteinIt->second.isBound())
        {
            complexPop.bindTo(secondProteinIt->second.getBindingSurface());
        }
    }
    
    // Apply dissociation if any occurs
    if (dissociatedAmount > 0 && complexIt != cell.m_molecules.end()) {
        // Remove from complex population
        complexIt->second.m_fNumber -= dissociatedAmount;
        assert(complexIt->second.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert protein level doesn't go below minimum
        
        // Return to free protein populations
        firstProteinIt->second.m_fNumber += dissociatedAmount;
        secondProteinIt->second.m_fNumber += dissociatedAmount;
    }
    
    return (boundAmount > 0 || dissociatedAmount > 0);
}
