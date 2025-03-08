#pragma once
#include <functional>
#include <memory>

class MRNA;
enum class CellCycleState;

class Organelle
{
public:
    // Base update function that takes cell cycle state and mRNA callback
    virtual void update(double dt, CellCycleState cellState, 
                       std::function<void(std::shared_ptr<MRNA>)> addMRNA) = 0;
    
    virtual ~Organelle() {}
};

