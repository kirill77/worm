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

void GridCell::updateMRNAs(double dt)
{
    // Handle mRNA degradation and cleanup
    auto it = m_pMRNAs.begin();
    while (it != m_pMRNAs.end()) {
        (*it)->degrade(dt); // Handle mRNA degradation via half-life
        if ((*it)->getNumber() <= 0.01) { // Remove degraded mRNAs
            it = m_pMRNAs.erase(it);
        } else {
            ++it;
        }
    }
} 