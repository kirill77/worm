#pragma once

#include <memory>
#include <vector>

// Forward declarations
class GPUMesh;
struct GPUMeshNode;

// this object is in 'gpu' and not in 'helpers' project because GPUWorld
// needs to know the declaration of this object
struct IVisObject
{
    virtual GPUMeshNode updateAndGetMeshNode() = 0;
}; 