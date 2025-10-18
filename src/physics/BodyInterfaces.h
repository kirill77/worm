#pragma once

#include <memory>
#include "geometry/vectors/vector.h"

class EdgeMesh;

// Per-vertex dynamic state (velocity, force, mass)
// Position is stored in the mesh geometry
struct INodeView
{
    double3 m_vVelocity;
    double3 m_vForce;
    double m_fMass;
};

// Interface for soft body physics operating on edge-based meshes
// Provides access to mesh topology (via m_pMesh) and per-vertex dynamic state (via getVertex)
struct IFaceBody : public std::enable_shared_from_this<IFaceBody>
{
    virtual ~IFaceBody() = default;

    std::shared_ptr<EdgeMesh> m_pMesh;

    virtual const INodeView& getVertex(uint32_t index) const = 0;
    virtual INodeView& getVertex(uint32_t index) = 0;
};

