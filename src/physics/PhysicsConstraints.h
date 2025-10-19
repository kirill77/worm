#pragma once

#include <memory>
#include "BodyInterfaces.h"

// Constraint interface (XPBD-style)
// Each constraint is bound to specific bodies at construction time
struct IConstraint : public std::enable_shared_from_this<IConstraint>
{
    virtual ~IConstraint() = default;
    virtual void project(double dt) = 0;
};


