#pragma once

#include <vector>
#include <array>
#include "math/vector.h"
#include "molecules/GridCell.h"

struct Grid
{
    // find cell by position
    GridCell& findCell(const float3& position);
    const GridCell& findCell(const float3& position) const;
    uint32_t positionToIndex(const float3& position) const;

    // find neighbors
    std::vector<uint32_t> getNeighborIndices(size_t cellIndex) const;

    // direct cell access
    size_t size() const { return m_grid.size(); }
    GridCell& operator[](size_t index) { return m_grid[index]; }
    const GridCell& operator[](size_t index) const { return m_grid[index]; }

private:
    static constexpr uint32_t GRID_RES = 3;  // 3x3x3 grid
    std::array<GridCell, GRID_RES* GRID_RES* GRID_RES> m_grid;
}; 