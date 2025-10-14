#pragma once

#include <memory>
#include "biology/simulation/TimeContext.h"
#include "biology/simulation/PhysicsCore.h"

class Cell;
class BVHMesh;
class PhysicsCore;

class CellSim
{
private:
    std::shared_ptr<Cell> m_pCell;
    std::shared_ptr<PhysicsCore> m_pPhysicsCore;

public:
    // Constructor
    CellSim(std::shared_ptr<Cell> pCell);
    
    // Simulation step method
    virtual void update(const TimeContext& time);
    
    // Accessor for the underlying cell
    std::shared_ptr<Cell> getCell() const { return m_pCell; }
};

