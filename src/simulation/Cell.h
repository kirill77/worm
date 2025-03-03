#pragma once

#include <vector>
#include <memory>

class Medium;

class Cell
{
private:
    std::vector<std::shared_ptr<class Organelle>> m_pOrganelles;
    std::shared_ptr<Medium> m_pMedium;

public:
    Cell(std::shared_ptr<Medium> pMedium);
    void simulateStep(double dt);
};

