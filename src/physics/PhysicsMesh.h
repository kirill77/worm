#pragma once

#include <memory>
#include <vector>
#include "geometry/vectors/vector.h"
#include "geometry/mesh/TriangleMesh.h"

// Per-vertex dynamic state (velocity, force, mass)
// Position is stored in the mesh geometry
struct PhysVertex
{
    double3 m_vVelocity;
    double3 m_vForce;
    double m_fMass;
};

// Physics mesh combining edge-based geometry with per-vertex dynamic state
class PhysicsMesh : public std::enable_shared_from_this<PhysicsMesh>
{
public:
    std::shared_ptr<TriangleMesh> m_pMesh;

    explicit PhysicsMesh(std::shared_ptr<TriangleMesh> mesh);

    const PhysVertex& getVertex(uint32_t index) const { return m_nodeData[index]; }
    PhysVertex& getVertex(uint32_t index) { return m_nodeData[index]; }

private:
    std::vector<PhysVertex> m_nodeData;
};

