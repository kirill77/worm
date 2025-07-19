#pragma once

#include <vector>
#include <memory>

class Organism
{
protected:
    std::vector<std::shared_ptr<class CellSim>> m_pCellSims;

public:
    virtual void simulateStep(double dt);

    const std::vector<std::shared_ptr<class CellSim>>& getCellSims() const
    {
        return m_pCellSims;
    }
};

