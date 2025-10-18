#pragma once

#include <memory>
#include "BodyInterfaces.h"

// Single-body constraint interface (XPBD-style)
struct IConstraint
{
    virtual ~IConstraint() = default;
    virtual void project(IFaceBody& body, double dt) = 0;
};

// Two-body constraint interface for interactions
struct ITwoBodyConstraint
{
    virtual ~ITwoBodyConstraint() = default;
    virtual void project(IFaceBody& a, IFaceBody& b, double dt) = 0;
};


