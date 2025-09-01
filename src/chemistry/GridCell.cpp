#include "pch.h"
#include "GridCell.h"
#include "Molecule.h"
#include "MoleculeWiki.h"
#include <cmath>

GridCell::GridCell()
{
    // Default constructor
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
    // Handle tRNA charging - attempt to charge uncharged tRNAs
    for (auto& pTRNA : m_pTRNAs) {
        pTRNA->charge(dt);
    }
    
    // Note: tRNAs don't degrade like mRNAs, so we don't remove them
    // They get recycled after being used in translation
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