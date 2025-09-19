#include "ComplexFormationInteraction.h"
#include "ResourceDistributor.h"
#include <algorithm>
#include <cmath>

ComplexFormationInteraction::ComplexFormationInteraction(
    const Molecule& firstProtein, 
    const Molecule& secondProtein, 
    const Parameters& params)
    : MoleculeInteraction(Mechanism::BINDING, 0.2)  // Lower ATP cost for binding
    , m_firstProtein(firstProtein)
    , m_secondProtein(secondProtein)
    , m_bindingRate(params.bindingRate)
    , m_dissociationRate(params.dissociationRate)
    , m_saturationConstant(params.saturationConstant)
    , m_complexId(params.complexId)
{
}

bool ComplexFormationInteraction::apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const
{
    double firstProteinAmount = resDistributor.getAvailableResource(m_firstProtein);
    double secondProteinAmount = resDistributor.getAvailableResource(m_secondProtein);
    
    // Calculate binding using mass action kinetics
    double bindingPotential = m_bindingRate * firstProteinAmount * secondProteinAmount / 
                             (m_saturationConstant + firstProteinAmount + secondProteinAmount);
    
    // The amount that can actually bind is limited by the lesser of the two proteins
    double boundAmount = std::min(bindingPotential * dt, std::min(firstProteinAmount, secondProteinAmount));
    
    // Binding requires ATP
    double requiredATP = boundAmount * m_atpCost;
    
    // Also check for dissociation of existing complexes (species-aware)
    Species species = m_firstProtein.getSpecies();
    // Both participants should share species in our loader; keep a defensive check
    assert(species == m_secondProtein.getSpecies());
    Molecule complexKey(m_complexId, ChemicalType::PROTEIN, species);
    auto complexIt = cell.m_molecules.find(complexKey);
    double complexAmount = (complexIt != cell.m_molecules.end()) ? complexIt->second.m_fNumber : 0.0;
    
    // Calculate dissociation (simpler first-order kinetics)
    double dissociatedAmount = complexAmount * m_dissociationRate * dt;
    
    // If we're in a dry run, just report resource requirements and return
    if (resDistributor.isDryRun()) {
        if (boundAmount > 0) {
            // Register our requirements with the resource distributor
            resDistributor.notifyResourceWanted(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE), requiredATP);
            resDistributor.notifyResourceWanted(m_firstProtein, boundAmount);
            resDistributor.notifyResourceWanted(m_secondProtein, boundAmount);
            return true; // We're reporting resource needs
        }
        // Dissociation doesn't consume resources, but still return true if it occurs
        return dissociatedAmount > 0;
    }

    auto firstProteinIt = cell.m_molecules.find(m_firstProtein);
    auto secondProteinIt = cell.m_molecules.find(m_secondProtein);

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
        auto& complexPop = cell.getOrCreateMolPop(complexKey);
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
