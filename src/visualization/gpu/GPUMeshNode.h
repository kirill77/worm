#pragma once

#include <memory>
#include <vector>
#include <DirectXMath.h>
#include "geometry/vectors/affine.h"
#include "geometry/vectors/box.h"

// Forward declarations
class GPUMesh;

/**
 * Represents a hierarchical node in a scene graph that can contain:
 * - Multiple meshes at this transform level
 * - Child nodes with their own relative transforms
 * This enables complex hierarchical objects and efficient instancing.
 */
struct GPUMeshNode
{
public:
    GPUMeshNode(const affine3& transform = affine3::identity());
    
    // Transform accessors
    const affine3& getTransform() const { return m_transform; }
    void setTransform(const affine3& transform) { m_transform = transform; }
    
    // Mesh management
    void addMesh(std::shared_ptr<GPUMesh> mesh);
    const std::vector<std::shared_ptr<GPUMesh>>& getMeshes() const { return m_meshes; }
    void clearMeshes() { m_meshes.clear(); }
    
    // Child node management
    void addChild(const GPUMeshNode& child);
    void addChild(GPUMeshNode&& child);
    const std::vector<GPUMeshNode>& getChildren() const { return m_children; }
    void clearChildren() { m_children.clear(); }
    
    // Convenience methods
    DirectX::XMMATRIX getWorldMatrix() const;
    DirectX::XMMATRIX getWorldMatrix(const affine3& parentTransform) const;
    box3 getWorldBoundingBox() const;
    box3 getWorldBoundingBox(const affine3& parentTransform) const;
    
    // Utility methods
    bool isEmpty() const { return m_meshes.empty() && m_children.empty(); }
    void clear() { clearMeshes(); clearChildren(); }
    
private:
    affine3 m_transform;  // Transform from node local space to parent space
    std::vector<std::shared_ptr<GPUMesh>> m_meshes;  // Meshes at this node level
    std::vector<GPUMeshNode> m_children;  // Child nodes with relative transforms
}; 