#pragma once
#include <memory>

class Cell;
class Medium;
enum class CellCycleState;

class Organelle
{
public:
    // Base update function that takes cell and medium
    virtual void update(double dt, Cell& cell, Medium& medium) = 0;
    
    virtual ~Organelle() {}
};

