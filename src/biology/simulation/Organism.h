#pragma once

#include <vector>
#include <memory>

class Organism
{
protected:
    std::vector<std::shared_ptr<class CellSim>> m_pCells;

public:
    virtual void simulateStep(double dt);

    const std::vector<std::shared_ptr<class CellSim>>& getCells() const
    {
        return m_pCells;
    }
};

