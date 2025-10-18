#pragma once

#include <memory>
#include "BodyInterfaces.h"

// Single-body constraint interface (XPBD-style)
// Each constraint is bound to specific bodies at construction time
struct IConstraint
{
    virtual ~IConstraint() = default;
    virtual void project(double dt) = 0;
};

// Two-body constraint interface for interactions
// Each constraint is bound to specific bodies at construction time
struct ITwoBodyConstraint
{
    virtual ~ITwoBodyConstraint() = default;
    virtual void project(double dt) = 0;
};


