#include "pch.h"
#include "Grid.h"
#include <cassert>

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