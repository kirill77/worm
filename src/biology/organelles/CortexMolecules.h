#pragma once

#include <unordered_map>
#include "geometry/mesh/MeshLocation.h"
#include "chemistry/molecules/Molecule.h"

// Cortex-bound molecular container (molecules bound at a surface location)
struct CortexMolecules : public MeshLocation
{
    std::unordered_map<Molecule, Population> m_bsMolecules;
};
