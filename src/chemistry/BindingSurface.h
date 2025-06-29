#pragma once

#include <memory>
#include <unordered_map>
#include <string>
#include "Molecule.h"
#include "geometry/vectors/vector.h"
#include "ProteinWiki.h"

/**
 * Abstract base class representing any cellular structure that can bind proteins.
 * This includes membranes, organelles, and other cellular structures.
 */
class BindingSurface : public std::enable_shared_from_this<BindingSurface>
{
protected:
    // Binding surface type
    StringDict::ID m_surfaceType;

public:
    /**
     * Constructor that initializes the binding surface
     * 
     * @param surfaceType Type of binding surface
     * @param bindingSiteDensity Amount of binding sites per surface unit (default: 1000)
     */
    BindingSurface(StringDict::ID surfaceType = StringDict::ID::eUNKNOWN)
        : m_surfaceType(surfaceType)
    {}
}; 