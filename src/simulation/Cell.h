#pragma once

#include <vector>
#include <memory>

class Cell
{
private:
    std::vector<std::unique_ptr<class Organelle>> m_pOrganelles;

public:
    Cell();
    void simulateStep(double dt);
};

