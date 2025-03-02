#include "pch.h"
#include "Worm.h"
#include "Cell.h"

Worm::Worm()
{
    auto pCell = std::make_shared<Cell>();
    m_pCells.push_back(pCell);
}
