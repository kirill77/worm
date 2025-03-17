#pragma once

#include <memory>
#include <unordered_map>
#include <string>
#include "Protein.h"
#include "math/vector.h"
#include "ProteinWiki.h"

/**
 * Abstract base class representing any cellular structure that can bind proteins.
 * This includes membranes, organelles, and other cellular structures.
 */
class ProteinBindingSurface : public std::enable_shared_from_this<ProteinBindingSurface>
{
protected:
    // Binding surface type
    ProteinWiki::BindingSurface m_surfaceType;

public:
    /**
     * Constructor that initializes the binding surface
     * 
     * @param surfaceType Type of binding surface
     * @param bindingSiteDensity Amount of binding sites per surface unit (default: 1000)
     */
    ProteinBindingSurface(ProteinWiki::BindingSurface surfaceType = ProteinWiki::BindingSurface::eUNKNOWN)
        : m_surfaceType(surfaceType)
    {}
}; 