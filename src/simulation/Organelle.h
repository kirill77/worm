#pragma once
#include <memory>

class Medium;
enum class CellCycleState;

class Organelle
{
public:
    // Base update function that takes cell cycle state and medium
    virtual void update(double dt, CellCycleState cellState, Medium& extMedium) = 0;
    
    virtual ~Organelle() {}
};

