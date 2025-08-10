#pragma once

#include <memory>
#include <vector>
#include "GPUMeshNode.h"

// Forward declarations
class GPUMesh;

// this object is in 'gpu' and not in 'helpers' project because GPUWorld
// needs to know the declaration of this object
struct IVisObject
{
public:
    // Non-virtual functions for caching mesh node
    GPUMeshNode updateMeshNode()
    {
        m_cachedMeshNode = updateAndGetMeshNode();
        return m_cachedMeshNode;
    }
    
    GPUMeshNode getMeshNode()
    {
        return m_cachedMeshNode;
    }
    
    // Virtual destructor for proper cleanup
    virtual ~IVisObject() = default;

protected:
    virtual GPUMeshNode updateAndGetMeshNode() = 0;

private:
    GPUMeshNode m_cachedMeshNode;
}; 