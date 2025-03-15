#pragma once

#include <memory>
#include "math/vector.h"

/**
 * Abstract base class representing any cellular structure that can bind proteins.
 * This includes membranes, organelles, and other cellular structures.
 */
class ProteinBindingSurface
{
protected:
    // Surface properties
    double m_fSurfaceArea;      // Surface area in square micrometers
    double m_fBindingCapacity;  // Maximum protein binding capacity per square micrometer

public:
    /**
     * Constructor that initializes the binding surface
     * 
     * @param fSurfaceArea Surface area in square micrometers
     * @param fBindingCapacity Maximum protein binding capacity per square micrometer
     */
    ProteinBindingSurface(double fSurfaceArea = 1.0, double fBindingCapacity = 1e6)
        : m_fSurfaceArea(fSurfaceArea)
        , m_fBindingCapacity(fBindingCapacity)
    {}
    
    /**
     * Virtual destructor for proper cleanup in derived classes
     */
    virtual ~ProteinBindingSurface() = default;
}; 