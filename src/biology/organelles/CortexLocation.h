#pragma once

#include <unordered_map>
#include <cstdint>
#include "geometry/vectors/vector.h"
#include "chemistry/molecules/Molecule.h"

struct Y_TuRC;

// Generic geometric address on a triangulated surface (e.g., cortex)
struct CortexLocation
{
    uint32_t m_triangleIndex = 0;
    float3 m_barycentric = float3(0,0,0);
    float3 m_normalized = float3(0,0,0); // normalized coordinates in cell medium space [-1,1]
};

// Cortex-bound molecular container (molecules bound at a surface location)
struct CortexMolecules : public CortexLocation
{
    std::unordered_map<Molecule, Population> m_bsMolecules;
};
