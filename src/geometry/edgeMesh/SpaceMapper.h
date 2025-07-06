#pragma once

#include <memory>
#include "geometry/vectors/vector.h"

struct EdgeMesh;

class SpaceMapper
{
    SpaceMapper(std::shared_ptr<EdgeMesh> pEdgeMesh);

    void update();

    // map from normalized space (-1, 1) to world space
    float3 normalizedToWorld(float3 vNormalized);

    // map from world space to normalized space (-1, 1)
    float3 worldToNormalized(float3 vWorld);

private:

};

