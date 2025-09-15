#include "PhosphorylationInteraction.h"
#include "chemistry/molecules/MoleculeWiki.h"
#include "ResourceDistributor.h"
#include <algorithm>
#include <cmath>

PhosphorylationInteraction::PhosphorylationInteraction(
    StringDict::ID kinaseId, 
    StringDict::ID targetId, 
    StringDict::ID phosphorylatedId,
    const Parameters& params)
    : MoleculeInteraction(Mechanism::PHOSPHORYLATION, 0.5)  // Standard ATP cost for phosphorylation
    , m_kinaseId(kinaseId)
    , m_targetId(targetId)
    , m_phosphorylatedId(phosphorylatedId)
    , m_removalRate(params.removalRate)
    , m_saturationConstant(params.saturationConstant)
{
}

bool PhosphorylationInteraction::apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const
{
    double kinaseAmount = resDistributor.getAvailableResource(Molecule(m_kinaseId, ChemicalType::PROTEIN));
    double targetAmount = resDistributor.getAvailableResource(Molecule(m_targetId, ChemicalType::PROTEIN));
    
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
            resDistributor.notifyResourceWanted(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE), requiredATP);
            resDistributor.notifyResourceWanted(Molecule(m_targetId, ChemicalType::PROTEIN), phosphorylatedAmount);
            return true; // We're reporting resource needs
        }
        return false;
    }

    // Apply the effect if any phosphorylation occurs
    if (phosphorylatedAmount > 0) {
        // Update ATP consumption
        auto& atpPop = cell.getOrCreateMolPop(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE));
        atpPop.m_fNumber -= requiredATP;
        assert(atpPop.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert ATP doesn't go below minimum
        
        // Remove proteins from unphosphorylated population
        auto targetIt = cell.m_molecules.find(Molecule(m_targetId, ChemicalType::PROTEIN));
        targetIt->second.m_fNumber -= phosphorylatedAmount;
        assert(targetIt->second.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert protein level doesn't go below minimum
        
        // Add to phosphorylated population
        auto& phosphorylatedPop = cell.getOrCreateMolPop(Molecule(m_phosphorylatedId, ChemicalType::PROTEIN));
        phosphorylatedPop.m_fNumber += phosphorylatedAmount;
        
        return true;
    }
    
    return false;
}
