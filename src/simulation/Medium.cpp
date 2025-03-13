#include "pch.h"
#include "Medium.h"
#include <random>
#include <algorithm>
#include <cassert>
#include "ProteinAntagonism.h"

ProteinPopulation& Medium::GridCell::getOrCreateProtein(const std::string& sProteinName)
{
    auto it = m_proteins.find(sProteinName);
    if (it != m_proteins.end()) {
        return it->second;
    }
    
    // Create new population with zero initial amount
    return m_proteins.emplace(sProteinName, ProteinPopulation(sProteinName, 0.0)).first->second;
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
    
    // Update each cell
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        bool isCortex = isCortexCell(i);
        if (!isCortex) continue;  // Skip non-cortical cells

        // Process each protein antagonism
        for (const auto& antagonism : m_proteinAntagonisms)
        {
            // Get protein amounts
            const auto& antagonist = antagonism.getAntagonist();
            const auto& target = antagonism.getTarget();
            
            double antagonistAmount = 0.0;
            double targetAmount = 0.0;
            
            // Get antagonist amount
            auto antagonistIt = m_grid[i].m_proteins.find(antagonist);
            if (antagonistIt != m_grid[i].m_proteins.end()) {
                antagonistAmount = antagonistIt->second.m_fNumber;
            }
            
            // Get target amount
            auto targetIt = m_grid[i].m_proteins.find(target);
            if (targetIt != m_grid[i].m_proteins.end()) {
                targetAmount = targetIt->second.m_fNumber;
            }
            
            // Calculate and apply antagonistic effects
            if (antagonistAmount > 0 && targetAmount > 0)
            {
                // Get available ATP if needed for phosphorylation
                double availableATP = 0.0;
                if (antagonism.getMechanism() == ProteinAntagonism::Mechanism::PHOSPHORYLATION) {
                    availableATP = m_grid[i].m_fAtp;
                }
                
                auto& targetPop = newGrid[i].getOrCreateProtein(target);
                
                // Calculate removal based on mechanism and parameters
                double removedAmount = antagonism.calculateRemoval(targetAmount, antagonistAmount, dt, availableATP);
                
                // Consume ATP if using phosphorylation
                if (antagonism.getMechanism() == ProteinAntagonism::Mechanism::PHOSPHORYLATION && removedAmount > 0) {
                    double atpUsed = removedAmount * antagonism.getATPCost();
                    newGrid[i].m_fAtp = std::max(0.0, newGrid[i].m_fAtp - atpUsed);
                }
                
                // Apply removal
                targetPop.m_fNumber = std::max(0.0, targetPop.m_fNumber - removedAmount);
                
                // Recovery of target
                double recoveredAmount = antagonism.calculateRecovery(removedAmount, dt);
                
                // Add recovered proteins to nearby non-cortical cells
                auto neighbors = getNeighborIndices(i);
                for (size_t neighborIdx : neighbors)
                {
                    if (!isCortexCell(neighborIdx))
                    {
                        auto& neighborPop = newGrid[neighborIdx].getOrCreateProtein(target);
                        neighborPop.m_fNumber += recoveredAmount / neighbors.size();
                    }
                }
            }
        }
        
        // Handle cortical binding/unbinding - this part remains similar to the original
        for (auto& [proteinName, pop] : m_grid[i].m_proteins)
        {
            if (!isCortex) // For non-cortical cells, this should not execute due to the skip above
            {
                // Skip, already filtered out
            }
            else // For cortical cells, handle proteins that can unbind
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
    }
    
    // Handle proteins in non-cortical cells that can bind to cortex
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        bool isCortex = isCortexCell(i);
        if (isCortex) continue; // Skip cortical cells, already processed
        
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
    // Initialize protein antagonisms
    
    // PKC-3 (kinase) phosphorylates posterior PARs
    ProteinAntagonism::Parameters pkc3ToParParams{
        0.9,                                          // High removal rate (strong kinase)
        0.07,                                         // Recovery rate
        550.0,                                        // Low saturation for stronger effect
        ProteinAntagonism::Mechanism::PHOSPHORYLATION, // Phosphorylation mechanism
        0.5                                           // ATP cost
    };
    
    // PAR-1 (kinase) phosphorylates PAR-3
    ProteinAntagonism::Parameters par1ToPar3Params{
        0.7,                                          // Medium-high removal rate
        0.06,                                         // Lower recovery rate
        650.0,                                        // Medium saturation constant
        ProteinAntagonism::Mechanism::PHOSPHORYLATION, // Phosphorylation mechanism
        0.4                                           // ATP cost
    };
    
    // PAR-1 weakly affects PAR-6 (indirect)
    ProteinAntagonism::Parameters par1ToPar6Params{
        0.3,                                          // Lower removal rate
        0.1,                                          // Higher recovery rate
        900.0,                                        // High saturation constant
        ProteinAntagonism::Mechanism::RECRUITMENT,     // Indirect mechanism
        0.0                                           // No ATP cost
    };
    
    // PAR-2 affects anterior proteins through cortical exclusion
    ProteinAntagonism::Parameters par2ToPar3Params{
        0.5,                                          // Medium removal rate 
        0.09,                                         // Medium recovery rate
        750.0,                                        // Medium saturation constant
        ProteinAntagonism::Mechanism::CORTICAL_EXCLUSION, // Cortical competition
        0.0                                           // No ATP cost
    };
    
    ProteinAntagonism::Parameters par2ToPar6Params{
        0.35,                                         // Medium-low removal rate
        0.09,                                         // Recovery rate
        800.0,                                        // Medium-high saturation constant
        ProteinAntagonism::Mechanism::CORTICAL_EXCLUSION, // Cortical competition
        0.0                                           // No ATP cost
    };
    
    ProteinAntagonism::Parameters par2ToPkc3Params{
        0.3,                                          // Lower removal rate
        0.1,                                          // Higher recovery rate
        850.0,                                        // Higher saturation constant
        ProteinAntagonism::Mechanism::CORTICAL_EXCLUSION, // Cortical competition
        0.0                                           // No ATP cost
    };
    
    // Add antagonistic relationships
    m_proteinAntagonisms.emplace_back("PKC-3", "PAR-2", pkc3ToParParams);
    m_proteinAntagonisms.emplace_back("PKC-3", "PAR-1", pkc3ToParParams);
    
    m_proteinAntagonisms.emplace_back("PAR-1", "PAR-3", par1ToPar3Params);
    m_proteinAntagonisms.emplace_back("PAR-1", "PAR-6", par1ToPar6Params);
    
    m_proteinAntagonisms.emplace_back("PAR-2", "PAR-3", par2ToPar3Params);
    m_proteinAntagonisms.emplace_back("PAR-2", "PAR-6", par2ToPar6Params);
    m_proteinAntagonisms.emplace_back("PAR-2", "PKC-3", par2ToPkc3Params);
}
