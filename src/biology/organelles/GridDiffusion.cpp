#include "pch.h"
#include "GridDiffusion.h"
#include "Grid.h"
#include <cassert>

GridDiffusion::GridDiffusion()
{
}

double GridDiffusion::computeDiffusionAmount(double moleculeCount, size_t numNeighbors, double dt) const
{
    return moleculeCount * DIFFUSION_RATE * dt / numNeighbors;
}

void GridDiffusion::updateDiffusion(Grid& grid, double dt)
{
    // Create an array of source populations that will participate in diffusion
    std::vector<uint32_t> nSourcePopsPerCell;
    std::vector<MPopulation*> pSourcePops;
    nSourcePopsPerCell.resize(grid.size(), 0);
    for (uint32_t uPass = 0; uPass < 2; ++uPass)
    {
        uint32_t nTotalSourcePops = 0;
        for (uint32_t uSourceCell = 0; uSourceCell < grid.size(); ++uSourceCell)
        {
            auto& sourcePops = grid[uSourceCell].m_molecules;
            for (auto itSourcePop = sourcePops.begin(); itSourcePop != sourcePops.end(); )
            {
                MPopulation& sourcePop = itSourcePop->second;
                if (sourcePop.m_population.m_fNumber == 0)
                {
                    itSourcePop = sourcePops.erase(itSourcePop);
                    continue;
                }
                ++itSourcePop;
                // if population is bound to a surface - it doesn't diffuse
                if (sourcePop.isBound())
                    continue;

                if (uPass == 0)
                {
                    ++nSourcePopsPerCell[uSourceCell];
                }
                else if (uPass == 1)
                {
                    pSourcePops[nTotalSourcePops] = &sourcePop;
                }
                ++nTotalSourcePops;
            }
        }

        if (uPass == 0)
        {
            pSourcePops.resize(nTotalSourcePops, nullptr);
        }
    }

    // Execute three passes: first pass allocates memory for diffusion amounts, second
    // pass calculation diffusion amount per destination population, and the third pass
    // deposits the diffused amount
    std::vector<double> diffusionAmounts;
    for (uint32_t uPass = 0; uPass < 3; ++uPass)
    {
        uint32_t uDiffusionIndex = 0;
        uint32_t uPopIndex = 0;

        for (uint32_t uSourceCell = 0; uSourceCell < grid.size(); ++uSourceCell)
        {
            auto vecNeighbors = grid.getNeighborIndices(uSourceCell);
            uint32_t uCellStartIndex = uPopIndex;

            if (uPass == 0)
            {
                uDiffusionIndex += nSourcePopsPerCell[uSourceCell] * (uint32_t)vecNeighbors.size();
            }
            else
            {
                // For each source population in this cell
                for (uint32_t uCellPopIndex = 0; uCellPopIndex < nSourcePopsPerCell[uSourceCell]; ++uCellPopIndex)
                {
                    MPopulation* pSourcePop = pSourcePops[uCellStartIndex + uCellPopIndex];
                    if (uPass == 1)
                    {
                        // Second pass: compute and store diffusion amounts
                        double fDiffusionAmount = computeDiffusionAmount(pSourcePop->m_population.m_fNumber, vecNeighbors.size(), dt);

                        // Store diffusion amounts for each neighbor
                        for (uint32_t uN = 0; uN < vecNeighbors.size(); ++uN)
                        {
                            diffusionAmounts[uDiffusionIndex + uN] = fDiffusionAmount;
                        }
                    }
                    else if (uPass == 2)
                    {
                        // Third pass: apply diffusion amounts
                        // Apply diffusion to each neighbor
                        for (uint32_t uN = 0; uN < vecNeighbors.size(); ++uN)
                        {
                            auto& destMPop = grid[vecNeighbors[uN]].getOrCreateMolecule(pSourcePop->getName());
                            pSourcePop->m_population.m_fNumber -= diffusionAmounts[uDiffusionIndex + uN];
                            destMPop.m_population.m_fNumber += diffusionAmounts[uDiffusionIndex + uN];
                        }
                    }
                    uDiffusionIndex += (uint32_t)vecNeighbors.size();
                }
            }
            uPopIndex += nSourcePopsPerCell[uSourceCell];
        }

        // After first pass, pre-allocate memory
        if (uPass == 0)
        {
            diffusionAmounts.resize(uDiffusionIndex);
        }
    }
} 