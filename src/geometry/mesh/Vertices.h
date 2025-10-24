#pragma once

#include <vector>
#include <cstdint>
#include "geometry/vectors/vector.h"
#include "geometry/vectors/box.h"
#include "Identifiable.h"

// Base mesh class containing only vertices
// Can represent point clouds or serve as base for 1D/2D/3D primitives
class Vertices : public Identifiable
{
public:
    static const uint32_t INVALID_INDEX = UINT32_MAX;
    
    struct Vertex {
        float3 position;
        Vertex(const float3& pos) : position(pos) {}
    };

    // Constructors
    Vertices();
    virtual ~Vertices() = default;
    
    // Vertex operations
    uint32_t addVertex(const float3& position);
    float3 getVertexPosition(uint32_t index) const;
    void setVertexPosition(uint32_t index, const float3& position);
    uint32_t getVertexCount() const;
    void removeLastVertex();
    
    // Clear mesh data
    virtual void clear();
    
    // Version tracking
    uint64_t getVersion() const { return m_version; }

    // Bounding box (cached based on version)
    box3 getBox() const;

protected:
    // Helper for derived classes to notify of changes
    void incrementVersion() { ++m_version; }

private:
    std::vector<Vertex> m_vertices;
    uint64_t m_version = 0;
    
    // Cached bounding box
    mutable box3 m_cachedBox = box3::empty();
    mutable uint64_t m_cachedBoxVersion = UINT64_MAX; // Invalid version initially
};

