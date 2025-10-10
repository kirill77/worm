#pragma once

#include <unordered_map>
#include <cstdint>
#include <cassert>
#include <limits>
#include <cmath>
#include "geometry/vectors/vector.h"
#include "chemistry/molecules/Molecule.h"

struct Y_TuRC;

// Generic geometric address on a triangulated surface (e.g., cortex)
struct CortexLocation
{
    uint32_t m_triangleIndex = 0;

private:
    float3 m_barycentric = float3(0,0,0);
    float3 m_normalized = float3(0,0,0); // normalized coordinates in cell medium space [-1,1]

    static float invalidSentinel() { return std::numeric_limits<float>::quiet_NaN(); }

public:
    const float3& getBarycentric() const { return m_barycentric; }
    void setBarycentric(const float3& v)
    {
        m_barycentric = v;
        // Updating barycentric invalidates cached normalized coordinate (debug-only sentinel)
#ifndef NDEBUG
        m_normalized.x = invalidSentinel();
#endif
    }

    const float3& getNormalized() const
    {
        assert(!std::isnan(m_normalized.x) && "Attempted to read invalid normalized coordinate; must be recomputed");
        return m_normalized;
    }
    void setNormalized(const float3& v)
    {
        assert(!std::isnan(v.x) && "Setting invalid normalized coordinate (NaN sentinel) is not allowed");
        m_normalized = v;
    }
};

// Cortex-bound molecular container (molecules bound at a surface location)
struct CortexMolecules : public CortexLocation
{
    std::unordered_map<Molecule, Population> m_bsMolecules;
};
