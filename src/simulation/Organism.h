#pragma once

#include <vector>
#include <memory>

class Organism
{
    std::vector<std::shared_ptr<class Cell>> m_pCells;

public:
    virtual void update(double dt) = 0;  // dt: time step
};

