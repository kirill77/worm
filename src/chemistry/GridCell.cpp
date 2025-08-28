#include "pch.h"
#include "GridCell.h"

GridCell::GridCell()
{
    // Default constructor
}

MPopulation& GridCell::getOrCreateMolecule(const std::string& sMoleculeName)
{
    auto it = m_molecules.find(sMoleculeName);
    if (it != m_molecules.end()) {
        return it->second;
    }
    
    // Create new population with zero initial amount
    return m_molecules.emplace(sMoleculeName, MPopulation(sMoleculeName, 0.0)).first->second;
}

MRNA& GridCell::getOrCreateMRNA(const std::string& sName)
{
    auto it = m_pMRNAs.find(sName);
    if (it != m_pMRNAs.end()) {
        return it->second;
    }
    
    // Create new mRNA with default parameters (will need to be set properly)
    return m_pMRNAs.emplace(sName, MRNA(sName, 0.0, 2.0, 1.0)).first->second;
}

void GridCell::updateMRNAs(double dt)
{
    // Handle mRNA degradation and cleanup
    auto it = m_pMRNAs.begin();
    while (it != m_pMRNAs.end()) {
        it->second.degrade(dt); // Handle mRNA degradation via half-life
        if (it->second.getNumber() <= 0.01) { // Remove degraded mRNAs
            it = m_pMRNAs.erase(it);
        } else {
            ++it;
        }
    }
}

void GridCell::updateTRNAs(double dt)
{
    // Handle tRNA charging - attempt to charge uncharged tRNAs
    for (auto& pTRNA : m_pTRNAs) {
        pTRNA->charge(dt);
    }
    
    // Note: tRNAs don't degrade like mRNAs, so we don't remove them
    // They get recycled after being used in translation
} 