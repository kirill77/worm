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
    // Structure to hold molecule identity and population reference for diffusion
    struct DiffusionEntry {
        const Molecule* molecule;
        Population* population;
        DiffusionEntry(const Molecule* mol, Population* pop) : molecule(mol), population(pop) {}
    };

    // Create arrays for efficient diffusion processing
    std::vector<uint32_t> nSourcePopsPerCell;
    std::vector<DiffusionEntry> sourcePops;
    nSourcePopsPerCell.resize(grid.size(), 0);
    
    for (uint32_t uPass = 0; uPass < 2; ++uPass)
    {
        uint32_t nTotalSourcePops = 0;
        for (uint32_t uSourceCell = 0; uSourceCell < grid.size(); ++uSourceCell)
        {
            auto& cellMolecules = grid[uSourceCell].m_molecules;
            for (auto itMolecule = cellMolecules.begin(); itMolecule != cellMolecules.end(); )
            {
                Population& population = itMolecule->second;
                if (population.m_fNumber == 0)
                {
                    itMolecule = cellMolecules.erase(itMolecule);
                    continue;
                }
                
                // if population is bound to a surface - it doesn't diffuse
                if (population.isBound()) {
                    ++itMolecule;
                    continue;
                }

                if (uPass == 0)
                {
                    ++nSourcePopsPerCell[uSourceCell];
                }
                else if (uPass == 1)
                {
                    sourcePops.emplace_back(&itMolecule->first, &population);
                }
                ++nTotalSourcePops;
                
                ++itMolecule;  // Only increment at the very end, after we've used the iterator
            }
        }

        if (uPass == 0)
        {
            sourcePops.reserve(nTotalSourcePops);
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
                    const DiffusionEntry& sourceEntry = sourcePops[uCellStartIndex + uCellPopIndex];
                    if (uPass == 1)
                    {
                        // Second pass: compute and store diffusion amounts
                        double fDiffusionAmount = computeDiffusionAmount(sourceEntry.population->m_fNumber, vecNeighbors.size(), dt);

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
                            auto& destPop = grid[vecNeighbors[uN]].getOrCreateMolPop(sourceEntry.molecule->getName());
                            sourceEntry.population->m_fNumber -= diffusionAmounts[uDiffusionIndex + uN];
                            destPop.m_fNumber += diffusionAmounts[uDiffusionIndex + uN];
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