#pragma once

#include <memory>
#include "BodyInterfaces.h"

// Single-body constraint interface (XPBD-style)
struct IConstraint
{
    virtual ~IConstraint() = default;
    virtual void project(IBody& body, double dt) = 0;
};

// Two-body constraint interface for interactions
struct ITwoBodyConstraint
{
    virtual ~ITwoBodyConstraint() = default;
    virtual void project(IBody& a, IBody& b, double dt) = 0;
};


