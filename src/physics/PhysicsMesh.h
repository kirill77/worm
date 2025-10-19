#pragma once

#include <memory>
#include <vector>
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

// Physics mesh combining edge-based geometry with per-vertex dynamic state
class PhysicsMesh : public std::enable_shared_from_this<PhysicsMesh>
{
public:
    std::shared_ptr<EdgeMesh> m_pMesh;

    explicit PhysicsMesh(std::shared_ptr<EdgeMesh> mesh);

    const INodeView& getVertex(uint32_t index) const { return m_nodeData[index]; }
    INodeView& getVertex(uint32_t index) { return m_nodeData[index]; }

private:
    std::vector<INodeView> m_nodeData;
};

