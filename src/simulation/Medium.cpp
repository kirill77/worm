#include "pch.h"
#include "Medium.h"
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

bool Medium::isCortexCell(size_t index) const
{
    // Convert to 3D coordinates
    size_t uZ = index % GRID_RES;
    size_t uY = (index / GRID_RES) % GRID_RES;
    size_t uX = index / (GRID_RES * GRID_RES);
    
    // Check if any coordinate is on the boundary
    return uX == 0 || uX == GRID_RES - 1 ||
           uY == 0 || uY == GRID_RES - 1 ||
           uZ == 0 || uZ == GRID_RES - 1;
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
    auto& gridCell = findCell(position);
    auto& proteinPop = gridCell.getOrCreateProtein(protein.m_sName);
    proteinPop.m_fNumber += protein.m_fNumber;
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

void Medium::updatePARDynamics(double fDt)
{
    // Create temporary grid for updated numbers
    auto gridNew = m_grid;
    
    // Constants for protein dynamics
    static constexpr double CORTEX_BINDING_RATE = 0.15;    // Rate at which proteins bind to cortex from cytoplasm
    static constexpr double CORTEX_UNBINDING_RATE = 0.12;  // Rate at which proteins unbind from cortex
    
    // Get all protein interactions 
    const auto& vecInteractions = ProteinWiki::GetProteinInteractions();
    
    // First, apply direct protein interactions
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        // Track ATP consumption
        double fAtpConsumed = 0.0;
        
        // Apply each interaction to the cell
        for (const auto& pInteraction : vecInteractions)
        {
            // Apply the interaction and track ATP consumption
            pInteraction->apply(gridNew[i], fDt, fAtpConsumed);
        }
        
        // Update ATP (already done inside apply() method, this is just for clarity)
        gridNew[i].m_fAtp = std::max(0.0, gridNew[i].m_fAtp);
    }
    
    // Then, handle neighbor effects for each cell
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        // Get references to neighbor cells
        auto vecNeighborIndices = getNeighborIndices(i);
        std::vector<std::reference_wrapper<GridCell>> vecNeighborCells;
        
        for (size_t uNeighborIdx : vecNeighborIndices)
        {
            vecNeighborCells.push_back(std::ref(gridNew[uNeighborIdx]));
        }
    }
    
    // Handle cortical binding/unbinding - this part remains similar to the original
    for (size_t i = 0; i < m_grid.size(); ++i)
    {
        bool bIsCortex = isCortexCell(i);
        
        if (bIsCortex)
        {
            // For cortical cells, handle proteins that can unbind
            for (auto& [sProteinName, proteinPop] : m_grid[i].m_proteins)
            {
                // Proteins can unbind from cortex to cytoplasm
                auto vecNeighbors = getNeighborIndices(i);
                std::vector<size_t> vecNonCortexNeighbors;
                
                // Find non-cortical neighbors
                for (size_t uNeighborIdx : vecNeighbors)
                {
                    if (!isCortexCell(uNeighborIdx))
                    {
                        vecNonCortexNeighbors.push_back(uNeighborIdx);
                    }
                }
                
                if (!vecNonCortexNeighbors.empty())
                {
                    // Calculate unbinding based on protein type
                    double fUnbindingRate = CORTEX_UNBINDING_RATE;
                    
                    // PAR-2 has more stability at the cortex
                    if (sProteinName == "PAR-2") {
                        fUnbindingRate *= 0.8; // 20% lower unbinding
                    }
                    
                    // Anterior complex proteins have slight stability
                    if (sProteinName == "PAR-3" || sProteinName == "PAR-6" || sProteinName == "PKC-3") {
                        fUnbindingRate *= 0.9; // 10% lower unbinding
                    }
                    
                    double fUnbindingAmount = proteinPop.m_fNumber * fUnbindingRate * fDt;
                    auto& proteinPopNew = gridNew[i].getOrCreateProtein(sProteinName);
                    proteinPopNew.m_fNumber -= fUnbindingAmount;
                    
                    // Distribute unbound proteins to non-cortical neighbors
                    double fAmountPerSite = fUnbindingAmount / vecNonCortexNeighbors.size();
                    for (size_t uNeighborIdx : vecNonCortexNeighbors)
                    {
                        auto& proteinPopNeighbor = gridNew[uNeighborIdx].getOrCreateProtein(sProteinName);
                        proteinPopNeighbor.m_fNumber += fAmountPerSite;
                    }
                }
            }
        }
        else
        {
            // For non-cortical cells, handle proteins that can bind to cortex
            auto vecNeighbors = getNeighborIndices(i);
            
            // For each protein population in the cell
            for (auto& [sProteinName, proteinPop] : m_grid[i].m_proteins)
            {
                std::vector<size_t> vecCortexNeighbors;
                
                // Find available cortex binding sites with protein-specific preferences
                for (size_t uNeighborIdx : vecNeighbors)
                {
                    if (isCortexCell(uNeighborIdx))
                    {
                        float3 vecPos = indexToPosition(uNeighborIdx);
                        
                        // For anterior PARs, prefer anterior cortex
                        if ((sProteinName == "PAR-3" || sProteinName == "PAR-6" || sProteinName == "PKC-3") &&
                            vecPos.y > 0)
                        {
                            vecCortexNeighbors.push_back(uNeighborIdx);
                        }
                        // For posterior PARs, prefer posterior cortex
                        else if ((sProteinName == "PAR-1" || sProteinName == "PAR-2") &&
                                vecPos.y < 0)
                        {
                            vecCortexNeighbors.push_back(uNeighborIdx);
                        }
                    }
                }
                
                if (!vecCortexNeighbors.empty())
                {
                    // Calculate binding amount with protein-specific rates
                    double fBindingRate = CORTEX_BINDING_RATE;
                    if (sProteinName == "PAR-2") fBindingRate *= 1.2; // 20% higher for PAR-2
                    if (sProteinName == "PAR-1") fBindingRate *= 1.1; // 10% higher for PAR-1
                    
                    double fBindingAmount = proteinPop.m_fNumber * fBindingRate * fDt;
                    auto& proteinPopNew = gridNew[i].getOrCreateProtein(sProteinName);
                    proteinPopNew.m_fNumber -= fBindingAmount;
                    
                    // Distribute to available cortex sites
                    double fAmountPerSite = fBindingAmount / vecCortexNeighbors.size();
                    for (size_t uNeighborIdx : vecCortexNeighbors)
                    {
                        auto& proteinPopNeighbor = gridNew[uNeighborIdx].getOrCreateProtein(sProteinName);
                        proteinPopNeighbor.m_fNumber += fAmountPerSite;
                    }
                }
            }
        }
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
    
    // Update PAR protein dynamics
    updatePARDynamics(fDt);
    
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
