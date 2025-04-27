#pragma once

#include <vector>
#include <array>
#include "math/vector.h"
#include "GridCell.h"

struct Grid
{
    // Helper functions
    GridCell& findCell(const float3& position);
    const GridCell& findCell(const float3& position) const;
    std::vector<uint32_t> getNeighborIndices(size_t cellIndex) const;

    // Add iterator support
    auto begin() { return m_grid.begin(); }
    auto end() { return m_grid.end(); }
    auto begin() const { return m_grid.begin(); }
    auto end() const { return m_grid.end(); }

    // Add size() function
    size_t size() const { return m_grid.size(); }

    // Add operator[] for direct access
    GridCell& operator[](size_t index) { return m_grid[index]; }
    const GridCell& operator[](size_t index) const { return m_grid[index]; }

private:
    // Convert between grid indices and 3D coordinates
    uint32_t positionToIndex(const float3& position) const;

    static constexpr uint32_t GRID_RES = 3;  // 3x3x3 grid
    std::array<GridCell, GRID_RES* GRID_RES* GRID_RES> m_grid;
}; 