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

    std::vector<double> diffusionAmounts(grid.size());

    // Execute three passes with the same code structure
    for (uint32_t pass = 0, uDiffusionIndex = 0; pass < 3; ++pass)
    {
        for (uint32_t i = 0; i < grid.size(); ++i)
        {
            auto vecNeighbors = grid.getNeighborIndices(i);

            // Process all proteins in the cell
            for (const auto& [proteinName, protein] : grid[i].m_proteins)
            {
                // Skip if protein is bound to a surface
                if (protein.isBound())
                    continue;

                if (pass == 1)
                {
                    // Second pass: compute and store diffusion amounts
                    double fDiffusionAmount = computeDiffusionAmount(protein.m_fNumber, vecNeighbors.size(), dt);
                    
                    // Store diffusion amounts for each neighbor
                    for (uint32_t uN = 0; uN < vecNeighbors.size(); ++uN)
                    {
                        diffusionAmounts[uDiffusionIndex + uN] = fDiffusionAmount;
                    }
                }
                else if (pass == 2)
                {
                    // Third pass: apply diffusion amounts
                    auto& mSource = grid[i].getOrCreateProtein(proteinName);
                    
                    // Store diffusion amounts for each neighbor
                    for (uint32_t uN = 0; uN < vecNeighbors.size(); ++uN)
                    {
                        auto& mDest = grid[vecNeighbors[uN]].getOrCreateProtein(proteinName);
                        mSource.m_fNumber -= diffusionAmounts[uDiffusionIndex + uN];
                        mDest.m_fNumber += diffusionAmounts[uDiffusionIndex + uN];
                    }
                }

                uDiffusionIndex += vecNeighbors.size();
            }

            // After first pass, pre-allocate memory
            if (pass == 0)
            {
                diffusionAmounts.resize(uDiffusionIndex);
            }
        }
    }

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