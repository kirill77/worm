#pragma once

#include <vector>
#include <memory>
#include "Medium.h"

enum class CellCycleState
{
    INTERPHASE,
    PROPHASE,
    METAPHASE,
    ANAPHASE,
    TELOPHASE,
    CYTOKINESIS
};

class Cell
{
private:
    std::vector<std::shared_ptr<class Organelle>> m_pOrganelles;
    std::shared_ptr<Medium> m_pMedium;
    CellCycleState m_cellCycleState;

public:
    Cell(std::shared_ptr<Medium> pMedium);
    void update(double dt);
    CellCycleState getCellCycleState() const { return m_cellCycleState; }
    std::shared_ptr<Medium> getMedium() const { return m_pMedium; }
};

