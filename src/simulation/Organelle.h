#pragma once
class Organelle
{
public:
    virtual void update(double dt) = 0;  // dt: time step
    virtual ~Organelle() {}
};

