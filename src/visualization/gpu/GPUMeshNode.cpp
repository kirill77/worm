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


box3 GPUMeshNode::getBoundingBox() const
{
    box3 result = box3::empty();
    bool hasValidBounds = false;

    // Accumulate bounding boxes from all meshes at this node level
    for (const auto& pMesh : m_meshes)
    {
        if (!pMesh)
            continue;

        const box3& localBounds = pMesh->getBoundingBox();
        if (!localBounds.isempty())
        {
            // Transform mesh bounds to this node's coordinate system
            box3 transformedBounds = localBounds * m_transform;

            if (hasValidBounds)
            {
                result = result | transformedBounds;
            }
            else
            {
                result = transformedBounds;
                hasValidBounds = true;
            }
        }
    }

    // Recursively accumulate bounding boxes from all child nodes
    for (const auto& child : m_children)
    {
        box3 childBounds = child.getBoundingBox();
        if (!childBounds.isempty())
        {
            // Child bounds are already in child's local space, transform to this node's space
            box3 transformedChildBounds = childBounds * m_transform;

            if (hasValidBounds)
            {
                result = result | transformedChildBounds;
            }
            else
            {
                result = transformedChildBounds;
                hasValidBounds = true;
            }
        }
    }

    return result;
}