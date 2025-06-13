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