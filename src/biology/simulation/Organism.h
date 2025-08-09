#pragma once

#include <vector>
#include <memory>
#include "biology/simulation/TimeContext.h"

class Organism
{
protected:
    std::vector<std::shared_ptr<class CellSim>> m_pCellSims;

public:
    virtual void simulateStep(const TimeContext& time);

    const std::vector<std::shared_ptr<class CellSim>>& getCellSims() const
    {
        return m_pCellSims;
    }
};

