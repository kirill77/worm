#include "pch.h"
#include "DephosphorylationInteraction.h"
#include "MoleculeWiki.h"
#include "ResourceDistributor.h"
#include <algorithm>

DephosphorylationInteraction::DephosphorylationInteraction(
    const std::string& target,
    const Parameters& params)
    : ProteinInteraction(Mechanism::DEPHOSPHORYLATION, 0.1)  // Lower ATP cost for dephosphorylation
    , m_targetName(target)
    , m_phosphorylatedName(MoleculeWiki::GetPhosphorylatedName(target))
    , m_recoveryRate(params.recoveryRate)
{
}

bool DephosphorylationInteraction::apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const
{
    // Calculate recovery
    double phosphorylatedAmount = resDistributor.getAvailableResource(m_phosphorylatedName);
    double recoveredAmount = phosphorylatedAmount * m_recoveryRate * dt;
    
    if (recoveredAmount <= 0) {
        return false;
    }
    
    // Dephosphorylation requires a small amount of ATP
    double requiredATP = recoveredAmount * m_atpCost;
    
    // If we're in a dry run, just report resource requirements and return
    if (resDistributor.isDryRun()) {
        // Register our requirements with the resource distributor
        resDistributor.notifyResourceWanted("ATP", requiredATP);
        resDistributor.notifyResourceWanted(m_phosphorylatedName, recoveredAmount);
        return true; // We're reporting resource needs
    }

    // Remove from phosphorylated population
    auto phosphorylatedIt = cell.m_molecules.find(Molecule(m_phosphorylatedName, ChemicalType::PROTEIN));
    phosphorylatedIt->second.m_fNumber -= recoveredAmount;
    assert(phosphorylatedIt->second.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert protein level doesn't go below minimum
    
    // Add back to original unphosphorylated population
    auto& unphosphorylatedPop = cell.getOrCreateMolPop(Molecule(m_targetName, ChemicalType::PROTEIN));
    unphosphorylatedPop.m_fNumber += recoveredAmount;
    
    // Update ATP consumption
    auto& atpPop = cell.getOrCreateMolPop(Molecule("ATP", ChemicalType::NUCLEOTIDE));
    atpPop.m_fNumber -= requiredATP;
    assert(atpPop.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert ATP doesn't go below minimum
    
    return true;
} 