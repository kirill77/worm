#include "GridCell.h"
#include "chemistry/molecules/Molecule.h"
#include "chemistry/molecules/MoleculeWiki.h"
#include "chemistry/molecules/TRNA.h"
#include <cmath>

GridCell::GridCell()
{
    // Default constructor
    m_volumeMicroM3 = 0.0;
}

Population& GridCell::getOrCreateMolPop(const Molecule& molecule)
{
    auto it = m_molecules.find(molecule);
    if (it != m_molecules.end()) {
        return it->second;
    }
    
    // Create new population with zero initial amount
    return m_molecules.emplace(molecule, Population(0.0)).first->second;
}

void GridCell::updateMRNAs(double dt)
{
    // Handle mRNA degradation and cleanup
    auto it = m_molecules.begin();
    while (it != m_molecules.end()) {
        if (it->first.getType() == ChemicalType::MRNA) {
            // Get half-life from MoleculeWiki
            const auto& info = MoleculeWiki::getInfo(it->first);
            double halfLife = info.m_fHalfLife;
            if (halfLife > 0.0) {
                // Simple exponential decay model for mRNA degradation
                it->second.m_fNumber *= exp(-dt / halfLife);
            }
            if (it->second.m_fNumber <= 0.01) { // Remove degraded mRNAs
                it = m_molecules.erase(it);
                continue;
            }
        }
        ++it;
    }
}

void GridCell::updateTRNAs(double dt)
{
    // Handle tRNA charging - convert uncharged tRNAs to charged tRNAs based on charging rate
    // Instead of iterating through all molecules, just check the known uncharged tRNA IDs
    const auto& unchargedTRNAIds = TRNA::getUnchargedTRNAIds();
    
    for (StringDict::ID unchargedID : unchargedTRNAIds) {
        Molecule unchargedTRNA(unchargedID, ChemicalType::TRNA);
        auto it = m_molecules.find(unchargedTRNA);
        
        if (it != m_molecules.end() && it->second.m_fNumber > 0.0) {
            const auto& info = MoleculeWiki::getInfo(it->first);
            double chargingRate = info.m_fChargingRate;
            
            if (chargingRate > 0.0) {
                // Calculate how many get charged in this time step
                double chargeProbability = chargingRate * dt;
                double chargedAmount = it->second.m_fNumber * chargeProbability;
                
                if (chargedAmount > 0.01) { // Only if significant amount
                    // Create corresponding charged tRNA molecule
                    StringDict::ID chargedID = TRNA::getChargedVariant(unchargedID);
                    Molecule chargedTRNA(chargedID, ChemicalType::TRNA);
                    
                    // Transfer molecules from uncharged to charged
                    Population& chargedPop = getOrCreateMolPop(chargedTRNA);
                    chargedPop.m_fNumber += chargedAmount;
                    it->second.m_fNumber -= chargedAmount;
                    
                    // Remove if fully consumed
                    if (it->second.m_fNumber <= 0.01) {
                        m_molecules.erase(it);
                    }
                }
            }
        }
    }
}

bool GridCell::hasMRNAs() const
{
    for (const auto& mol : m_molecules) {
        if (mol.first.getType() == ChemicalType::MRNA) {
            return true;
        }
    }
    return false;
}