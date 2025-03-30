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

void Medium::updateProteinDiffusionGrid(double dt)
{
    // Original grid-based diffusion implementation
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

void Medium::updateProteinDiffusionPhysical(double dt)
{
    // Create temporary grid for updated state
    auto gridNew = m_grid;
    
    // Process a fixed number of random samples
    for (int sample = 0; sample < DIFFUSION_SAMPLES; ++sample)
    {
        // Generate a random position in the medium
        float3 sourcePos = generateRandomPosition();
        
        // Get the grid cell at this position
        size_t sourceIndex = positionToIndex(sourcePos);
        GridCell& sourceCell = m_grid[sourceIndex];
        
        // Skip if there are no proteins in this cell
        if (sourceCell.m_proteins.empty()) {
            continue;
        }
        
        // Randomly select a protein population from this cell
        // First, build a vector of protein names with non-zero amounts
        std::vector<std::string> availableProteins;
        for (const auto& [proteinName, population] : sourceCell.m_proteins) {
            // Only include unbound proteins with positive amounts
            if (population.m_fNumber > 0.001 && !population.isBound()) {
                availableProteins.push_back(proteinName);
            }
        }
        
        // Skip if no eligible proteins
        if (availableProteins.empty()) {
            continue;
        }
        
        // Select a random protein
        std::uniform_int_distribution<size_t> proteinDist(0, availableProteins.size() - 1);
        const std::string& selectedProtein = availableProteins[proteinDist(g_rng)];
        
        // Get the protein population
        const auto& proteinPop = sourceCell.m_proteins.at(selectedProtein);
        double proteinAmount = proteinPop.m_fNumber;
        
        // Generate random direction and distance
        float3 direction = generateRandomDirection();
        float distance = generateRandomDistance(dt);
        
        // Calculate the new position
        float3 newPos = sourcePos + direction * distance;
        
        // Clamp to boundaries [-1, 1]
        newPos.x = std::clamp(newPos.x, -1.0f, 1.0f);
        newPos.y = std::clamp(newPos.y, -1.0f, 1.0f);
        newPos.z = std::clamp(newPos.z, -1.0f, 1.0f);
        
        // Get the target cell index
        size_t targetIndex = positionToIndex(newPos);
        
        // Skip if same cell (no movement)
        if (targetIndex == sourceIndex) {
            continue;
        }
        
        // Determine amount to transfer
        // Use a fraction based on distance (longer distances transfer less)
        // This prevents too many proteins from moving too far at once
        double maxTransferFraction = 0.05; // Max 5% per transfer
        double transferFraction = maxTransferFraction / (1.0 + distance * 5.0);
        double transferAmount = proteinAmount * transferFraction;
        
        // Apply a minimum threshold to avoid very small transfers
        if (transferAmount < 0.001) {
            continue;
        }

        // Update protein amounts in source and target cells
        // Use the getOrCreateProtein method instead of direct [] access to avoid default constructor issues
        ProteinPopulation& sourceProtein = gridNew[sourceIndex].getOrCreateProtein(selectedProtein);
        sourceProtein.m_fNumber -= transferAmount;
        
        ProteinPopulation& targetProtein = gridNew[targetIndex].getOrCreateProtein(selectedProtein);
        targetProtein.m_fNumber += transferAmount;
    }
    
    m_grid = std::move(gridNew);
}

void Medium::updateProteinDiffusion(double dt)
{
    // Choose which diffusion method to use
    // Uncomment the desired method:
    
    // Original grid-based method (only adjacent cells)
    //updateProteinDiffusionGrid(dt);
    
    // New physically-based method (arbitrary distances)
    updateProteinDiffusionPhysical(dt);
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
        for (int i = 0; i < vecInteractions.size(); ++i)
        {
            m_resDistributor.notifyNewInteractionStarting(*vecInteractions[i]);
            vecInteractions[i]->apply(m_grid[uCell], fDt, m_resDistributor);
        }

        // now do the real run to distribute the resources
        m_resDistributor.notifyNewRealRun();
        for (int i = 0; i < vecInteractions.size(); ++i)
        {
            if (!m_resDistributor.notifyNewInteractionStarting(*vecInteractions[i]))
                continue;
            vecInteractions[i]->apply(m_grid[uCell], fDt, m_resDistributor);
        }

        // Ensure ATP doesn't go below zero
        m_grid[uCell].m_fAtp = std::max(0.0, m_grid[uCell].m_fAtp);
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
