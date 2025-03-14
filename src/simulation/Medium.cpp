#include "pch.h"
#include "Medium.h"
#include <random>
#include <algorithm>
#include <cassert>
#include "GridCell.h"

uint32_t Medium::positionToIndex(const float3& position) const
{
    uint32_t index = 0;
    for (int i = 0; i < 3; ++i)
    {
        assert(position[i] >= -1.0f && position[i] <= 1.0f);
        float normalized = (position[i] + 1.0f) / 2.0f;
        uint32_t gridPos = std::min(GRID_RES - 1, static_cast<uint32_t>(GRID_RES * normalized));
        index = (index * GRID_RES) + gridPos;
    }
    return index;
}

float3 Medium::indexToPosition(size_t index) const
{
    float3 pos;
    for (int i = 2; i >= 0; --i)
    {
        size_t gridPos = index % GRID_RES;
        index /= GRID_RES;
        pos[i] = (2.0f * gridPos / (GRID_RES - 1.0f)) - 1.0f;
    }
    return pos;
}

GridCell& Medium::findCell(const float3& position)
{
    return m_grid[positionToIndex(position)];
}

const GridCell& Medium::findCell(const float3& position) const
{
    return m_grid[positionToIndex(position)];
}

bool Medium::isCortexCell(size_t index) const
{
    // Convert to 3D coordinates
    size_t z = index % GRID_RES;
    size_t y = (index / GRID_RES) % GRID_RES;
    size_t x = index / (GRID_RES * GRID_RES);
    
    // Check if any coordinate is on the boundary
    return x == 0 || x == GRID_RES - 1 ||
           y == 0 || y == GRID_RES - 1 ||
           z == 0 || z == GRID_RES - 1;
}

std::vector<size_t> Medium::getNeighborIndices(size_t cellIndex) const
{
    std::vector<size_t> neighbors;
    size_t z = cellIndex % GRID_RES;
    size_t y = (cellIndex / GRID_RES) % GRID_RES;
    size_t x = cellIndex / (GRID_RES * GRID_RES);
    
    // Check all 6 face neighbors
    const int offsets[6][3] = {
        {-1, 0, 0}, {1, 0, 0},
        {0, -1, 0}, {0, 1, 0},
        {0, 0, -1}, {0, 0, 1}
    };
    
    for (const auto& offset : offsets)
    {
        int newX = static_cast<int>(x) + offset[0];
        int newY = static_cast<int>(y) + offset[1];
        int newZ = static_cast<int>(z) + offset[2];
        
        if (newX >= 0 && newX < GRID_RES &&
            newY >= 0 && newY < GRID_RES &&
            newZ >= 0 && newZ < GRID_RES)
        {
            size_t neighborIndex = newX * GRID_RES * GRID_RES + newY * GRID_RES + newZ;
            neighbors.push_back(neighborIndex);
        }
    }
    
    return neighbors;
}

void Medium::addProtein(const ProteinPopulation& protein, const float3& position)
{
    auto& cell = findCell(position);
    auto& pop = cell.getOrCreateProtein(protein.m_sName);
    pop.m_fNumber += protein.m_fNumber;
}

void Medium::addMRNA(std::shared_ptr<MRNA> mRNA, const float3& position)
{
    findCell(position).m_pMRNAs.push_back(mRNA);
}

double Medium::getProteinNumber(const std::string& proteinName, const float3& position) const
{
    const auto& cell = findCell(position);
    auto it = cell.m_proteins.find(proteinName);
    return (it != cell.m_proteins.end()) ? it->second.m_fNumber : 0.0;
}

void Medium::updateProteinDiffusion(double dt)
{
    // Create temporary grid for updated numbers
    auto newGrid = m_grid;
    
    // Update each cell
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        auto neighbors = getNeighborIndices(i);
        
        // For each protein population in the cell
        for (auto& [proteinName, pop] : m_grid[i].m_proteins)
        {
            // Calculate amount to diffuse
            double diffusionAmount = pop.m_fNumber * DIFFUSION_RATE * dt / neighbors.size();
            
            // Get reference to population in new grid for this cell
            auto& sourcePop = newGrid[i].getOrCreateProtein(proteinName);
            
            // Distribute to neighbors
            for (size_t neighborIdx : neighbors)
            {
                auto& neighborPop = newGrid[neighborIdx].getOrCreateProtein(proteinName);
                neighborPop.m_fNumber += diffusionAmount;
                sourcePop.m_fNumber -= diffusionAmount;
            }
        }
    }
    
    m_grid = std::move(newGrid);
}

void Medium::updatePARDynamics(double dt)
{
    // Create temporary grid for updated numbers
    auto newGrid = m_grid;
    
    // Constants for protein dynamics
    static constexpr double CORTEX_BINDING_RATE = 0.15;    // Rate at which proteins bind to cortex from cytoplasm
    static constexpr double CORTEX_UNBINDING_RATE = 0.12;  // Rate at which proteins unbind from cortex
    
    // Get all protein interactions 
    const auto& interactions = ProteinWiki::GetProteinInteractions();
    
    // First, apply direct protein interactions
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        bool isCortex = isCortexCell(i);
        
        // Skip non-cortical cells for certain interactions
        if (!isCortex) continue;
        
        // Track ATP consumption
        double atpConsumed = 0.0;
        
        // Apply each interaction to the cell
        for (const auto& interaction : interactions)
        {
            // Apply the interaction and track ATP consumption
            interaction->apply(newGrid[i], dt, atpConsumed);
        }
        
        // Update ATP (already done inside apply() method, this is just for clarity)
        newGrid[i].m_fAtp = std::max(0.0, newGrid[i].m_fAtp);
    }
    
    // Then, handle neighbor effects for each cell
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        // Get references to neighbor cells
        auto neighborIndices = getNeighborIndices(i);
        std::vector<std::reference_wrapper<GridCell>> neighborCells;
        
        for (size_t neighborIdx : neighborIndices)
        {
            neighborCells.push_back(std::ref(newGrid[neighborIdx]));
        }
        
        // Apply neighbor effects for each interaction
        for (const auto& interaction : interactions)
        {
            interaction->applyNeighborEffects(newGrid[i], neighborCells, dt);
        }
    }
    
    // Handle cortical binding/unbinding - this part remains similar to the original
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        bool isCortex = isCortexCell(i);
        
        if (isCortex)
        {
            // For cortical cells, handle proteins that can unbind
            for (auto& [proteinName, pop] : m_grid[i].m_proteins)
            {
                // Proteins can unbind from cortex to cytoplasm
                auto neighbors = getNeighborIndices(i);
                std::vector<size_t> nonCortexNeighbors;
                
                // Find non-cortical neighbors
                for (size_t neighborIdx : neighbors)
                {
                    if (!isCortexCell(neighborIdx))
                    {
                        nonCortexNeighbors.push_back(neighborIdx);
                    }
                }
                
                if (!nonCortexNeighbors.empty())
                {
                    // Calculate unbinding based on protein type
                    double unbindingRate = CORTEX_UNBINDING_RATE;
                    
                    // PAR-2 has more stability at the cortex
                    if (proteinName == "PAR-2") {
                        unbindingRate *= 0.8; // 20% lower unbinding
                    }
                    
                    // Anterior complex proteins have slight stability
                    if (proteinName == "PAR-3" || proteinName == "PAR-6" || proteinName == "PKC-3") {
                        unbindingRate *= 0.9; // 10% lower unbinding
                    }
                    
                    double unbindingAmount = pop.m_fNumber * unbindingRate * dt;
                    auto& newPop = newGrid[i].getOrCreateProtein(proteinName);
                    newPop.m_fNumber -= unbindingAmount;
                    
                    // Distribute unbound proteins to non-cortical neighbors
                    double amountPerSite = unbindingAmount / nonCortexNeighbors.size();
                    for (size_t neighborIdx : nonCortexNeighbors)
                    {
                        auto& neighborPop = newGrid[neighborIdx].getOrCreateProtein(proteinName);
                        neighborPop.m_fNumber += amountPerSite;
                    }
                }
            }
        }
        else
        {
            // For non-cortical cells, handle proteins that can bind to cortex
            auto neighbors = getNeighborIndices(i);
            
            // For each protein population in the cell
            for (auto& [proteinName, pop] : m_grid[i].m_proteins)
            {
                std::vector<size_t> cortexNeighbors;
                
                // Find available cortex binding sites with protein-specific preferences
                for (size_t neighborIdx : neighbors)
                {
                    if (isCortexCell(neighborIdx))
                    {
                        float3 pos = indexToPosition(neighborIdx);
                        
                        // For anterior PARs, prefer anterior cortex
                        if ((proteinName == "PAR-3" || proteinName == "PAR-6" || proteinName == "PKC-3") &&
                            pos.y > 0)
                        {
                            cortexNeighbors.push_back(neighborIdx);
                        }
                        // For posterior PARs, prefer posterior cortex
                        else if ((proteinName == "PAR-1" || proteinName == "PAR-2") &&
                                pos.y < 0)
                        {
                            cortexNeighbors.push_back(neighborIdx);
                        }
                    }
                }
                
                if (!cortexNeighbors.empty())
                {
                    // Calculate binding amount with protein-specific rates
                    double bindingRate = CORTEX_BINDING_RATE;
                    if (proteinName == "PAR-2") bindingRate *= 1.2; // 20% higher for PAR-2
                    if (proteinName == "PAR-1") bindingRate *= 1.1; // 10% higher for PAR-1
                    
                    double bindingAmount = pop.m_fNumber * bindingRate * dt;
                    auto& newPop = newGrid[i].getOrCreateProtein(proteinName);
                    newPop.m_fNumber -= bindingAmount;
                    
                    // Distribute to available cortex sites
                    double amountPerSite = bindingAmount / cortexNeighbors.size();
                    for (size_t neighborIdx : cortexNeighbors)
                    {
                        auto& neighborPop = newGrid[neighborIdx].getOrCreateProtein(proteinName);
                        neighborPop.m_fNumber += amountPerSite;
                    }
                }
            }
        }
    }
    
    m_grid = std::move(newGrid);
}

double Medium::getTotalProteinNumber(const std::string& proteinName) const
{
    double total = 0.0;
    for (const auto& cell : m_grid)
    {
        auto it = cell.m_proteins.find(proteinName);
        if (it != cell.m_proteins.end()) {
            total += it->second.m_fNumber;
        }
    }
    return total;
}

void Medium::update(double dt)
{
    // Update diffusion of proteins and ATP
    updateProteinDiffusion(dt);
    updateATPDiffusion(dt);
    
    // Update PAR protein dynamics
    updatePARDynamics(dt);
    
    // Update mRNA positions
    translateMRNAs(dt);
}

void Medium::translateMRNAs(double dt)
{
    // TODO: Implement translation of mRNAs into proteins
    // This will need to:
    // 1. Check for available tRNAs
    // 2. Create new proteins
    // 3. Add proteins to appropriate cytoplasmic regions
}

void Medium::addATP(double amount, const float3& position)
{
    auto& cell = findCell(position);
    cell.m_fAtp = std::min(cell.m_fAtp + amount, MAX_ATP_PER_CELL);
}

bool Medium::consumeATP(double amount, const float3& position)
{
    auto& cell = findCell(position);
    if (cell.m_fAtp >= amount)
    {
        cell.m_fAtp -= amount;
        return true;
    }
    return false;
}

double Medium::getAvailableATP(const float3& position) const
{
    const auto& cell = findCell(position);
    return cell.m_fAtp;
}

void Medium::updateATPDiffusion(double dt)
{
    // Create a temporary copy of ATP levels
    std::vector<double> newATPLevels(m_grid.size());
    
    // Calculate diffusion for each cell
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        double currentATP = m_grid[i].m_fAtp;
        auto neighbors = getNeighborIndices(i);
        
        // Calculate net ATP change due to diffusion
        double atpChange = 0.0;
        for (size_t neighborIdx : neighbors)
        {
            double neighborATP = m_grid[neighborIdx].m_fAtp;
            double diffusion = (neighborATP - currentATP) * ATP_DIFFUSION_RATE * dt;
            atpChange += diffusion;
        }
        
        // Store new ATP level
        newATPLevels[i] = std::min(MAX_ATP_PER_CELL, 
                                  std::max(0.0, currentATP + atpChange));
    }
    
    // Update ATP levels
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        m_grid[i].m_fAtp = newATPLevels[i];
    }
}

Medium::Medium()
{
    // No need to initialize protein antagonisms here anymore
    // They are now managed by ProteinWiki
}
