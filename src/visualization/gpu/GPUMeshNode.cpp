#include "GPUMeshNode.h"
#include "GPUMesh.h"

GPUMeshNode::GPUMeshNode(const affine3& transform)
    : m_transform(transform)
{
}

void GPUMeshNode::addMesh(std::shared_ptr<GPUMesh> mesh)
{
    if (mesh)
    {
        m_meshes.push_back(mesh);
    }
}

void GPUMeshNode::addChild(const GPUMeshNode& child)
{
    m_children.push_back(child);
}

void GPUMeshNode::addChild(GPUMeshNode&& child)
{
    m_children.push_back(std::move(child));
}

DirectX::XMMATRIX GPUMeshNode::getWorldMatrix() const
{
    return getWorldMatrix(affine3::identity());
}

DirectX::XMMATRIX GPUMeshNode::getWorldMatrix(const affine3& parentTransform) const
{
    // Combine parent transform with this node's transform
    affine3 worldTransform = parentTransform * m_transform;
    
    // Convert affine3 to DirectX::XMMATRIX
    const auto& m = worldTransform.m_linear;
    const auto& t = worldTransform.m_translation;
    
    return DirectX::XMMATRIX(
        m.m00, m.m01, m.m02, 0.0f,
        m.m10, m.m11, m.m12, 0.0f,
        m.m20, m.m21, m.m22, 0.0f,
        t.x,   t.y,   t.z,   1.0f
    );
}

box3 GPUMeshNode::getWorldBoundingBox() const
{
    return getWorldBoundingBox(affine3::identity());
}

box3 GPUMeshNode::getWorldBoundingBox(const affine3& parentTransform) const
{
    // Combine parent transform with this node's transform
    affine3 worldTransform = parentTransform * m_transform;
    
    box3 boundingBox = box3::empty();
    bool hasValidBounds = false;
    
    // Include bounding boxes of all meshes at this node level
    for (const auto& mesh : m_meshes)
    {
        if (mesh)
        {
            const box3& localBounds = mesh->getBoundingBox();
            if (!localBounds.isempty())
            {
                // Transform mesh bounds to world space
                box3 worldBounds = localBounds * worldTransform;
                
                if (hasValidBounds)
                {
                    boundingBox = boundingBox | worldBounds;
                }
                else
                {
                    boundingBox = worldBounds;
                    hasValidBounds = true;
                }
            }
        }
    }
    
    // Include bounding boxes of all child nodes
    for (const auto& child : m_children)
    {
        box3 childBounds = child.getWorldBoundingBox(worldTransform);
        if (!childBounds.isempty())
        {
            if (hasValidBounds)
            {
                boundingBox = boundingBox | childBounds;
            }
            else
            {
                boundingBox = childBounds;
                hasValidBounds = true;
            }
        }
    }
    
    return boundingBox;
} 