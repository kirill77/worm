#include "pch.h"
#include "Medium.h"
#include "ProteinWiki.h"
#include "ResourceAllocation.h"
#include <random>
#include <algorithm>
#include <cassert>
#include "GridCell.h"

uint32_t Medium::positionToIndex(const float3& position) const
{
    uint32_t uIndex = 0;
    for (int i = 0; i < 3; ++i)
    {
        assert(position[i] >= -1.0f && position[i] <= 1.0f);
        float fNormalized = (position[i] + 1.0f) / 2.0f;
        uint32_t uGridPos = std::min(GRID_RES - 1, static_cast<uint32_t>(GRID_RES * fNormalized));
        uIndex = (uIndex * GRID_RES) + uGridPos;
    }
    return uIndex;
}

float3 Medium::indexToPosition(size_t index) const
{
    float3 vecPos;
    for (int i = 2; i >= 0; --i)
    {
        size_t uGridPos = index % GRID_RES;
        index /= GRID_RES;
        vecPos[i] = (2.0f * uGridPos / (GRID_RES - 1.0f)) - 1.0f;
    }
    return vecPos;
}

GridCell& Medium::findCell(const float3& position)
{
    return m_grid[positionToIndex(position)];
}

const GridCell& Medium::findCell(const float3& position) const
{
    return m_grid[positionToIndex(position)];
}

std::vector<size_t> Medium::getNeighborIndices(size_t cellIndex) const
{
    std::vector<size_t> vecNeighbors;
    size_t uZ = cellIndex % GRID_RES;
    size_t uY = (cellIndex / GRID_RES) % GRID_RES;
    size_t uX = cellIndex / (GRID_RES * GRID_RES);
    
    // Check all 6 face neighbors
    const int aOffsets[6][3] = {
        {-1, 0, 0}, {1, 0, 0},
        {0, -1, 0}, {0, 1, 0},
        {0, 0, -1}, {0, 0, 1}
    };
    
    for (const auto& offset : aOffsets)
    {
        int iNewX = static_cast<int>(uX) + offset[0];
        int iNewY = static_cast<int>(uY) + offset[1];
        int iNewZ = static_cast<int>(uZ) + offset[2];
        
        if (iNewX >= 0 && iNewX < GRID_RES &&
            iNewY >= 0 && iNewY < GRID_RES &&
            iNewZ >= 0 && iNewZ < GRID_RES)
        {
            size_t uNeighborIndex = iNewX * GRID_RES * GRID_RES + iNewY * GRID_RES + iNewZ;
            vecNeighbors.push_back(uNeighborIndex);
        }
    }
    
    return vecNeighbors;
}

void Medium::addProtein(const ProteinPopulation& protein, const float3& position)
{
    GridCell& gridCell = findCell(position);
    ProteinPopulation& cellProtein = gridCell.getOrCreateProtein(protein.m_sName);

    cellProtein.bindTo(protein.getBindingSurface());
    cellProtein.m_fNumber += protein.m_fNumber;
}

void Medium::addMRNA(std::shared_ptr<MRNA> pMRNA, const float3& position)
{
    findCell(position).m_pMRNAs.push_back(pMRNA);
}

double Medium::getProteinNumber(const std::string& proteinName, const float3& position) const
{
    const auto& gridCell = findCell(position);
    auto itProtein = gridCell.m_proteins.find(proteinName);
    return (itProtein != gridCell.m_proteins.end()) ? itProtein->second.m_fNumber : 0.0;
}

void Medium::updateProteinDiffusion(double dt)
{
    // Create temporary grid for updated numbers
    auto gridNew = m_grid;
    
    // Update each cell
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        auto vecNeighbors = getNeighborIndices(i);
        
        // For each protein population in the cell
        for (auto& [sProteinName, proteinPop] : m_grid[i].m_proteins)
        {
            // proteins attached to surfaces don't participate in diffusion
            if (proteinPop.isBound())
                continue;

            // Calculate amount to diffuse
            double fDiffusionAmount = proteinPop.m_fNumber * DIFFUSION_RATE * dt / vecNeighbors.size();
            
            // Get reference to population in new grid for this cell
            auto& proteinPopSource = gridNew[i].getOrCreateProtein(sProteinName);
            
            // Distribute to neighbors
            for (size_t uNeighborIdx : vecNeighbors)
            {
                auto& proteinPopNeighbor = gridNew[uNeighborIdx].getOrCreateProtein(sProteinName);
                proteinPopNeighbor.m_fNumber += fDiffusionAmount;
                proteinPopSource.m_fNumber -= fDiffusionAmount;
            }
        }
    }
    
    m_grid = std::move(gridNew);
}

void Medium::updateProteinInteraction(double fDt)
{
    // Create temporary grid for updated numbers
    auto gridNew = m_grid;
    
    // Get all protein interactions 
    const auto& vecInteractions = ProteinWiki::GetProteinInteractions();
    
    // First, apply direct protein interactions
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        // at first make a dry run to figure out who needs which resources
        m_resDistributor.notifyNewDryRun(m_grid[i]);
        for (int i = 0; i < vecInteractions.size(); ++i)
        {
            vecInteractions[i]->apply(gridNew[i], fDt, m_resDistributor);
        }

        // now do the real run to distribute the resources
        m_resDistributor.notifyNewRealRun();
        for (int i = 0; i < vecInteractions.size(); ++i)
        {
            vecInteractions[i]->apply(gridNew[i], fDt, m_resDistributor);
        }

        // Ensure ATP doesn't go below zero
        gridNew[i].m_fAtp = std::max(0.0, gridNew[i].m_fAtp);
    }

    m_grid = std::move(gridNew);
}

double Medium::getTotalProteinNumber(const std::string& proteinName) const
{
    double fTotal = 0.0;
    for (const auto& gridCell : m_grid)
    {
        auto itProtein = gridCell.m_proteins.find(proteinName);
        if (itProtein != gridCell.m_proteins.end()) {
            fTotal += itProtein->second.m_fNumber;
        }
    }
    return fTotal;
}

void Medium::update(double fDt)
{
    // Update diffusion of proteins and ATP
    updateProteinDiffusion(fDt);
    updateATPDiffusion(fDt);
    
    // Interaction of proteins between each other
    updateProteinInteraction(fDt);
    
    // Update mRNA positions
    translateMRNAs(fDt);
}

void Medium::translateMRNAs(double fDt)
{
    // TODO: Implement translation of mRNAs into proteins
    // This will need to:
    // 1. Check for available tRNAs
    // 2. Create new proteins
    // 3. Add proteins to appropriate cytoplasmic regions
}

void Medium::addATP(double fAmount, const float3& position)
{
    auto& gridCell = findCell(position);
    gridCell.m_fAtp = std::min(gridCell.m_fAtp + fAmount, MAX_ATP_PER_CELL);
}

bool Medium::consumeATP(double fAmount, const float3& position)
{
    auto& gridCell = findCell(position);
    if (gridCell.m_fAtp >= fAmount)
    {
        gridCell.m_fAtp -= fAmount;
        return true;
    }
    return false;
}

double Medium::getAvailableATP(const float3& position) const
{
    const auto& gridCell = findCell(position);
    return gridCell.m_fAtp;
}

void Medium::updateATPDiffusion(double fDt)
{
    // Create a temporary copy of ATP levels
    std::vector<double> vecNewATPLevels(m_grid.size());
    
    // Calculate diffusion for each cell
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        double fCurrentATP = m_grid[i].m_fAtp;
        auto vecNeighbors = getNeighborIndices(i);
        
        // Calculate net ATP change due to diffusion
        double fAtpChange = 0.0;
        for (size_t uNeighborIdx : vecNeighbors)
        {
            double fNeighborATP = m_grid[uNeighborIdx].m_fAtp;
            double fDiffusion = (fNeighborATP - fCurrentATP) * ATP_DIFFUSION_RATE * fDt;
            fAtpChange += fDiffusion;
        }
        
        // Store new ATP level
        vecNewATPLevels[i] = std::min(MAX_ATP_PER_CELL, 
                                  std::max(0.0, fCurrentATP + fAtpChange));
    }
    
    // Update ATP levels
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        m_grid[i].m_fAtp = vecNewATPLevels[i];
    }
}

Medium::Medium()
{
    // No need to initialize protein antagonisms here anymore
    // They are now managed by ProteinWiki
}
