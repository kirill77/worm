#pragma once

#include <memory>

class Cell;

class CellSim
{
private:
    std::shared_ptr<Cell> m_pCell;

public:
    // Constructor
    CellSim(std::shared_ptr<Cell> pCell);
    
    // Simulation step method
    virtual void update(double dt);
    
    // Accessor for the underlying cell
    std::shared_ptr<Cell> getCell() const { return m_pCell; }
};

