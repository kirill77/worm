#include "pch.h"
#include "Medium.h"
#include "ProteinWiki.h"
#include "ResourceDistributor.h"
#include <random>
#include <algorithm>
#include <cassert>
#include "GridCell.h"

// Global random number generator for consistent randomness
static std::mt19937 g_rng(std::random_device{}());

uint32_t Grid::positionToIndex(const float3& position) const
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

float3 Grid::indexToPosition(size_t index) const
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

GridCell& Grid::findCell(const float3& position)
{
    return m_grid[positionToIndex(position)];
}

const GridCell& Grid::findCell(const float3& position) const
{
    return m_grid[positionToIndex(position)];
}

std::vector<uint32_t> Grid::getNeighborIndices(size_t cellIndex) const
{
    std::vector<uint32_t> vecNeighbors;
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
            uint32_t uNeighborIndex = static_cast<uint32_t>(iNewX * GRID_RES * GRID_RES + iNewY * GRID_RES + iNewZ);
            vecNeighbors.push_back(uNeighborIndex);
        }
    }
    
    return vecNeighbors;
}

void Medium::addProtein(const MPopulation& protein, const float3& position)
{
    GridCell& gridCell = m_grid.findCell(position);
    MPopulation& cellProtein = gridCell.getOrCreateProtein(protein.m_sName);

    cellProtein.bindTo(protein.getBindingSurface());
    cellProtein.m_fNumber += protein.m_fNumber;
}

void Medium::addMRNA(std::shared_ptr<MRNA> pMRNA, const float3& position)
{
    m_grid.findCell(position).m_pMRNAs.push_back(pMRNA);
}

double Medium::getProteinNumber(const std::string& proteinName, const float3& position) const
{
    const auto& gridCell = m_grid.findCell(position);
    auto itProtein = gridCell.m_proteins.find(proteinName);
    return (itProtein != gridCell.m_proteins.end()) ? itProtein->second.m_fNumber : 0.0;
}

float3 Medium::generateRandomPosition() const
{
    // Generate random coordinates in [-1,1] range
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    return float3(dist(g_rng), dist(g_rng), dist(g_rng));
}

float3 Medium::generateRandomDirection() const
{
    // Generate random direction vector on unit sphere
    std::normal_distribution<float> dist(0.0f, 1.0f);
    float3 dir(dist(g_rng), dist(g_rng), dist(g_rng));
    
    // Normalize to get unit vector
    float len = length(dir);
    if (len > 0.001f) {
        dir = dir * (1.0f / len);
    } else {
        // Fallback for unlikely case of very small vector
        dir = float3(1.0f, 0.0f, 0.0f);
    }
    
    return dir;
}

float Medium::generateRandomDistance(double dt) const
{
    // Use normal distribution for distance to match diffusion physics
    // Scale by diffusion rate and time step (sqrt of dt for correct diffusion scaling)
    float sigma = static_cast<float>(DIFFUSION_SIGMA * std::sqrt(dt * DIFFUSION_RATE));
    std::normal_distribution<float> dist(0.0f, sigma);
    
    // Take absolute value to ensure positive distance
    return std::abs(dist(g_rng));
}

void Medium::updateProteinDiffusion(double dt)
{
    // Original grid-based diffusion implementation
    // Create temporary grid for updated numbers
    auto gridNew = m_grid;

    // Update each cell
    for (uint32_t i = 0; i < m_grid.getNCells(); ++i)
    {
        auto vecNeighbors = m_grid.getNeighborIndices(i);

        // For each protein population in the cell
        for (auto& [sProteinName, proteinPop] : m_grid.getCell(i).m_proteins)
        {
            // proteins attached to surfaces don't participate in diffusion
            if (proteinPop.isBound())
                continue;

            // Calculate amount to diffuse
            double fDiffusionAmount = proteinPop.m_fNumber * DIFFUSION_RATE * dt / vecNeighbors.size();

            // Get reference to population in new grid for this cell
            auto& proteinPopSource = gridNew.getCell(i).getOrCreateProtein(sProteinName);

            // Distribute to neighbors
            for (uint32_t uNeighborIdx : vecNeighbors)
            {
                auto& proteinPopNeighbor = gridNew.getCell(uNeighborIdx).getOrCreateProtein(sProteinName);
                proteinPopNeighbor.m_fNumber += fDiffusionAmount;
                proteinPopSource.m_fNumber -= fDiffusionAmount;
            }
        }
    }

    m_grid = std::move(gridNew);
}

void Medium::updateProteinInteraction(double fDt)
{
    // Get all protein interactions 
    const auto& vecInteractions = ProteinWiki::GetProteinInteractions();
    
    // First, apply direct protein interactions
    for (size_t uCell = 0; uCell < m_grid.size(); ++uCell)
    {
        // at first make a dry run to figure out who needs which resources
        m_resDistributor.notifyNewDryRun(m_grid[uCell]);
        for (size_t i = 0; i < vecInteractions.size(); ++i)
        {
            m_resDistributor.notifyNewInteractionStarting(*vecInteractions[i]);
            vecInteractions[i]->apply(m_grid[uCell], fDt, m_resDistributor);
        }

        // now do the real run to distribute the resources
        m_resDistributor.notifyNewRealRun();
        for (size_t i = 0; i < vecInteractions.size(); ++i)
        {
            if (!m_resDistributor.notifyNewInteractionStarting(*vecInteractions[i]))
                continue;
            vecInteractions[i]->apply(m_grid[uCell], fDt, m_resDistributor);
        }

        // Ensure ATP doesn't go below zero
        m_grid[uCell].m_fAtp = std::max<double>(0.0, m_grid[uCell].m_fAtp);
    }
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
    auto& gridCell = m_grid.findCell(position);
    gridCell.m_fAtp = std::min<double>(gridCell.m_fAtp + fAmount, MAX_ATP_PER_CELL);
}

bool Medium::consumeATP(double fAmount, const float3& position)
{
    auto& gridCell = m_grid.findCell(position);
    if (gridCell.m_fAtp >= fAmount)
    {
        gridCell.m_fAtp -= fAmount;
        return true;
    }
    return false;
}

double Medium::getAvailableATP(const float3& position) const
{
    const auto& gridCell = m_grid.findCell(position);
    return gridCell.m_fAtp;
}

void Medium::updateATPDiffusion(double fDt)
{
    // Create temporary grid for updated numbers
    auto gridNew = m_grid;

    // Update each cell
    for (uint32_t i = 0; i < m_grid.getNCells(); ++i)
    {
        auto vecNeighbors = m_grid.getNeighborIndices(i);
        
        // Calculate amount to diffuse
        double fDiffusionAmount = m_grid[i].m_fAtp * ATP_DIFFUSION_RATE * fDt / vecNeighbors.size();
        
        // Distribute to neighbors
        for (uint32_t uNeighborIdx : vecNeighbors)
        {
            gridNew[uNeighborIdx].m_fAtp += fDiffusionAmount;
            gridNew[i].m_fAtp -= fDiffusionAmount;
        }
    }

    m_grid = std::move(gridNew);
}

Medium::Medium()
{
    // No need to initialize protein antagonisms here anymore
    // They are now managed by ProteinWiki
}
