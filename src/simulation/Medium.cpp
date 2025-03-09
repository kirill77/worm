#include "pch.h"
#include "Medium.h"
#include <random>
#include <algorithm>
#include <cassert>

ProteinPopulation& Medium::GridCell::findOrCreatePopulation(const ProteinPopulation& protein)
{
    auto it = std::find_if(m_proteins.begin(), m_proteins.end(),
        [&](const ProteinPopulation& pop) {
            return pop.m_sName == protein.m_sName;
        });
    
    if (it != m_proteins.end()) {
        it->m_fNumber += protein.m_fNumber;  // Add to existing population
        return *it;
    }
    
    m_proteins.push_back(protein);  // Create new population
    return m_proteins.back();
}

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

Medium::GridCell& Medium::findCell(const float3& position)
{
    return m_grid[positionToIndex(position)];
}

const Medium::GridCell& Medium::findCell(const float3& position) const
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
    cell.findOrCreatePopulation(protein);
}

void Medium::addMRNA(std::shared_ptr<MRNA> mRNA, const float3& position)
{
    findCell(position).m_pMRNAs.push_back(mRNA);
}

double Medium::getProteinNumber(const std::string& proteinName, const float3& position) const
{
    const auto& cell = findCell(position);
    auto it = std::find_if(cell.m_proteins.begin(), cell.m_proteins.end(),
        [&](const ProteinPopulation& pop) {
            return pop.m_sName == proteinName;
        });
    
    return (it != cell.m_proteins.end()) ? it->m_fNumber : 0.0;
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
        for (auto& pop : m_grid[i].m_proteins)
        {
            double diffusionAmount = pop.m_fNumber * DIFFUSION_RATE * dt / neighbors.size();
            
            // Distribute to neighbors
            for (size_t neighborIdx : neighbors)
            {
                auto& neighborPop = newGrid[neighborIdx].findOrCreatePopulation(pop);
                neighborPop.m_fNumber += diffusionAmount;
                newGrid[i].findOrCreatePopulation(pop).m_fNumber -= diffusionAmount;
            }
        }
    }
    
    m_grid = std::move(newGrid);
}

void Medium::updatePARDynamics(double dt)
{
    // Create temporary grid for updated numbers
    auto newGrid = m_grid;
    
    // Update each cell
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        bool isCortex = isCortexCell(i);
        
        for (auto& pop : m_grid[i].m_proteins)
        {
            if (!isPARProtein(pop)) continue;
            
            if (isCortex)
            {
                // Handle unbinding from cortex
                double unbindingAmount = pop.m_fNumber * CORTEX_UNBINDING_RATE * dt;
                auto& newPop = newGrid[i].findOrCreatePopulation(pop);
                newPop.m_fNumber -= unbindingAmount;
                
                // Move to neighboring internal cell
                auto neighbors = getNeighborIndices(i);
                for (size_t neighborIdx : neighbors)
                {
                    if (!isCortexCell(neighborIdx))
                    {
                        auto& neighborPop = newGrid[neighborIdx].findOrCreatePopulation(pop);
                        neighborPop.m_fNumber += unbindingAmount / neighbors.size();
                    }
                }
            }
            else
            {
                // Handle binding to cortex
                double bindingAmount = pop.m_fNumber * CORTEX_BINDING_RATE * dt;
                auto& newPop = newGrid[i].findOrCreatePopulation(pop);
                newPop.m_fNumber -= bindingAmount;
                
                // Move to neighboring cortex cells
                auto neighbors = getNeighborIndices(i);
                std::vector<size_t> cortexNeighbors;
                for (size_t neighborIdx : neighbors)
                {
                    if (isCortexCell(neighborIdx))
                    {
                        cortexNeighbors.push_back(neighborIdx);
                    }
                }
                
                if (!cortexNeighbors.empty())
                {
                    double amountPerNeighbor = bindingAmount / cortexNeighbors.size();
                    for (size_t neighborIdx : cortexNeighbors)
                    {
                        auto& neighborPop = newGrid[neighborIdx].findOrCreatePopulation(pop);
                        neighborPop.m_fNumber += amountPerNeighbor;
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
        for (const auto& pop : cell.m_proteins)
        {
            if (pop.m_sName == proteinName)
            {
                total += pop.m_fNumber;
            }
        }
    }
    return total;
}

bool Medium::isPARProtein(const Protein& protein) const
{
    const std::string& name = protein.m_sName;
    return name.substr(0, 4) == "PAR-" || name == "PKC-3";
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
