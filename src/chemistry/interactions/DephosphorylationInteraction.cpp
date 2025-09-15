#include "DephosphorylationInteraction.h"
#include "chemistry/molecules/MoleculeWiki.h"
#include "ResourceDistributor.h"
#include <algorithm>

DephosphorylationInteraction::DephosphorylationInteraction(
    StringDict::ID targetId,
    StringDict::ID phosphorylatedId,
    const Parameters& params)
    : MoleculeInteraction(Mechanism::DEPHOSPHORYLATION, 0.1)  // Lower ATP cost for dephosphorylation
    , m_targetId(targetId)
    , m_phosphorylatedId(phosphorylatedId)
    , m_recoveryRate(params.recoveryRate)
{
}

bool DephosphorylationInteraction::apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const
{
    // Calculate recovery
    double phosphorylatedAmount = resDistributor.getAvailableResource(Molecule(m_phosphorylatedId, ChemicalType::PROTEIN));
    double recoveredAmount = phosphorylatedAmount * m_recoveryRate * dt;
    
    if (recoveredAmount <= 0) {
        return false;
    }
    
    // Dephosphorylation requires a small amount of ATP
    double requiredATP = recoveredAmount * m_atpCost;
    
    // If we're in a dry run, just report resource requirements and return
    if (resDistributor.isDryRun()) {
        // Register our requirements with the resource distributor
        resDistributor.notifyResourceWanted(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE), requiredATP);
        resDistributor.notifyResourceWanted(Molecule(m_phosphorylatedId, ChemicalType::PROTEIN), recoveredAmount);
        return true; // We're reporting resource needs
    }

    // Remove from phosphorylated population
    auto phosphorylatedIt = cell.m_molecules.find(Molecule(m_phosphorylatedId, ChemicalType::PROTEIN));
    phosphorylatedIt->second.m_fNumber -= recoveredAmount;
    assert(phosphorylatedIt->second.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert protein level doesn't go below minimum
    
    // Add back to original unphosphorylated population
    auto& unphosphorylatedPop = cell.getOrCreateMolPop(Molecule(m_targetId, ChemicalType::PROTEIN));
    unphosphorylatedPop.m_fNumber += recoveredAmount;
    
    // Update ATP consumption
            auto& atpPop = cell.getOrCreateMolPop(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE));
    atpPop.m_fNumber -= requiredATP;
    assert(atpPop.m_fNumber >= GridCell::MIN_RESOURCE_LEVEL); // Assert ATP doesn't go below minimum
    
    return true;
} 