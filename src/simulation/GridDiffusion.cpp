#include "pch.h"
#include "GridDiffusion.h"
#include "Grid.h"
#include <cassert>

GridDiffusion::GridDiffusion(const Parameters& params)
    : m_params(params)
{
}

double GridDiffusion::computeDiffusionAmount(double moleculeCount, size_t numNeighbors, double dt) const
{
    return moleculeCount * m_params.diffusionRate * dt / numNeighbors;
}

void GridDiffusion::updateDiffusion(Grid& grid, const std::string& moleculeName, double dt)
{
    copyATPToProteins(grid);

    // Create temporary grid for updated numbers
    auto gridNew = grid;

    // Update each cell
    for (uint32_t i = 0; i < grid.size(); ++i)
    {
        auto vecNeighbors = grid.getNeighborIndices(i);

        // Find the protein population in the cell
        auto itProtein = grid[i].m_proteins.find(moleculeName);
        if (itProtein == grid[i].m_proteins.end())
            continue;

        // Skip if protein is bound to a surface
        if (itProtein->second.isBound())
            continue;

        // Calculate amount to diffuse
        double fDiffusionAmount = computeDiffusionAmount(itProtein->second.m_fNumber, vecNeighbors.size(), dt);

        // Get reference to population in new grid for this cell
        auto& proteinPopSource = gridNew[i].getOrCreateProtein(moleculeName);

        // Distribute to neighbors
        for (uint32_t uNeighborIdx : vecNeighbors)
        {
            auto& proteinPopNeighbor = gridNew[uNeighborIdx].getOrCreateProtein(moleculeName);
            proteinPopNeighbor.m_fNumber += fDiffusionAmount;
            proteinPopSource.m_fNumber -= fDiffusionAmount;
        }
    }

    grid = std::move(gridNew);

    copyATPFromProteins(grid);
}

void GridDiffusion::copyATPToProteins(Grid& grid) const
{
    const std::string ATP_PROTEIN_NAME = "ATP";
    
    // For each cell in the grid
    for (uint32_t i = 0; i < grid.size(); ++i)
    {
        // Get or create the ATP protein population
        auto& atpProtein = grid[i].getOrCreateProtein(ATP_PROTEIN_NAME);
        
        // Copy the ATP value from m_fAtp to the protein population
        atpProtein.m_fNumber = grid[i].m_fAtp;
    }
}

void GridDiffusion::copyATPFromProteins(Grid& grid) const
{
    const std::string ATP_PROTEIN_NAME = "ATP";
    
    // For each cell in the grid
    for (uint32_t i = 0; i < grid.size(); ++i)
    {
        // Find the ATP protein population
        auto itProtein = grid[i].m_proteins.find(ATP_PROTEIN_NAME);
        
        // If ATP protein exists, copy its value to m_fAtp
        if (itProtein != grid[i].m_proteins.end())
        {
            grid[i].m_fAtp = itProtein->second.m_fNumber;
        }
        else
        {
            // If no ATP protein exists, set ATP to zero
            grid[i].m_fAtp = 0.0;
        }
    }
} 