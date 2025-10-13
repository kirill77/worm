#pragma once

#include <memory>
#include "physics/tensionSphere/BodyInterfaces.h"
#include "geometry/mesh/edgeMesh.h"

// Adapter that exposes EdgeMesh + per-vertex state as IBody views
class SoftBodyMeshAdapter : public IBody
{
public:
    SoftBodyMeshAdapter(std::shared_ptr<EdgeMesh> mesh,
        std::vector<double3>& velocities,
        const std::vector<double>& edgeRestLengths)
        : m_mesh(std::move(mesh))
        , m_velocities(velocities)
        , m_edgeRestLengths(edgeRestLengths)
    {}

    INodeView& nodes() override { return m_nodeView; }
    const INodeView& nodes() const override { return m_nodeView; }
    IEdgeView& edges() override { return m_edgeView; }
    const IEdgeView& edges() const override { return m_edgeView; }
    IFaceView& faces() override { return m_faceView; }
    const IFaceView& faces() const override { return m_faceView; }

private:
    struct NodeViewImpl : public INodeView
    {
        SoftBodyMeshAdapter& parent;
        explicit NodeViewImpl(SoftBodyMeshAdapter& p) : parent(p) {}
        size_t size() const override { return parent.m_mesh->getVertexCount(); }
        double3 getPosition(uint32_t i) const override { return double3(parent.m_mesh->getVertexPosition(i)); }
        double3 getVelocity(uint32_t i) const override { return parent.m_velocities[i]; }
        double  getMass(uint32_t i) const override { (void)i; return 1.0; }
        void    setPosition(uint32_t i, const double3& p) override { parent.m_mesh->setVertexPosition(i, float3(p)); }
        void    setVelocity(uint32_t i, const double3& v) override { parent.m_velocities[i] = v; }
        void    addForce(uint32_t i, const double3& f) override { parent.m_forces[i] += f; }
    } m_nodeView { *this };

    struct EdgeViewImpl : public IEdgeView
    {
        SoftBodyMeshAdapter& parent;
        explicit EdgeViewImpl(SoftBodyMeshAdapter& p) : parent(p) {}
        size_t edgeCount() const override { return parent.m_mesh->getEdgeCount(); }
        std::pair<uint32_t,uint32_t> edge(uint32_t e) const override { return parent.m_mesh->getEdge(e); }
        double restLength(uint32_t e) const override { return parent.m_edgeRestLengths[e]; }
    } m_edgeView { *this };

    struct FaceViewImpl : public IFaceView
    {
        SoftBodyMeshAdapter& parent;
        explicit FaceViewImpl(SoftBodyMeshAdapter& p) : parent(p) {}
        size_t faceCount() const override { return parent.m_mesh->getTriangleCount(); }
        uint3  face(uint32_t f) const override { return parent.m_mesh->getTriangleVertices(f); }
    } m_faceView { *this };

public:
    // Force buffer management
    void resizeForces()
    {
        m_forces.assign(m_mesh->getVertexCount(), double3(0,0,0));
    }
    std::vector<double3>& forces() { return m_forces; }

private:
    std::shared_ptr<EdgeMesh> m_mesh;
    std::vector<double3>& m_velocities;
    const std::vector<double>& m_edgeRestLengths;
    std::vector<double3> m_forces;
};


