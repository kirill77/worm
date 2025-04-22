#include "pch.h"
#include "GridCell.h"

GridCell::GridCell() : m_fAtp(0.0) 
{
    // Default constructor initializes ATP to zero
}

MPopulation& GridCell::getOrCreateProtein(const std::string& sProteinName)
{
    auto it = m_proteins.find(sProteinName);
    if (it != m_proteins.end()) {
        return it->second;
    }
    
    // Create new population with zero initial amount
    return m_proteins.emplace(sProteinName, MPopulation(sProteinName, 0.0)).first->second;
} 