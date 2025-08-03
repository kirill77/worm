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
