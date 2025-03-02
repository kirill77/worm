#pragma once

#include <vector>
#include <memory>

class Organism
{
protected:
    std::vector<std::shared_ptr<class Cell>> m_pCells;

public:
    virtual void simulateStep(double dt);
};

