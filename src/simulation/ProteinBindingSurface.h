#pragma once

#include <memory>
#include <unordered_map>
#include <string>
#include "Protein.h"
#include "math/vector.h"

/**
 * Abstract base class representing any cellular structure that can bind proteins.
 * This includes membranes, organelles, and other cellular structures.
 */
class ProteinBindingSurface : public std::enable_shared_from_this<ProteinBindingSurface>
{
protected:
    // Map of bound proteins and their quantities
    std::unordered_map<std::string, std::shared_ptr<ProteinPopulation>> m_boundProteins;

public:
    /**
     * Constructor that initializes the binding surface
     * 
     * @param fSurfaceArea Surface area in square micrometers
     * @param fBindingCapacity Maximum protein binding capacity per square micrometer
     */
    ProteinBindingSurface() {}
}; 