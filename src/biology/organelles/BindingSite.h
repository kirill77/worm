#pragma once

#include <unordered_map>
#include <cstdint>
#include "geometry/vectors/vector.h"
#include "chemistry/molecules/Molecule.h"

struct BindingSite
{
	uint32_t m_triangleIndex;
	float3 m_barycentric;
	float3 m_normalized; // normalized coordinates in cell medium space [-1,1]
	std::unordered_map<Molecule, Population> m_bsMolecules;
};


