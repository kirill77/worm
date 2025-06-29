#pragma once

#include <vector>
#include <memory>

struct Organism
{
protected:
    std::vector<std::shared_ptr<class Cell>> m_pCells;

public:
    virtual void simulateStep(double dt);

    const std::vector<std::shared_ptr<class Cell>>& getCells() const
    {
        return m_pCells;
    }
};

