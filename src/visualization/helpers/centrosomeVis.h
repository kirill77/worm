#pragma once

#include "visualization/gpu/IVisObject.h"
#include <memory>

struct Centrosome;
class GPUMesh;
struct GPUQueue;

struct CentrosomeVis : public IVisObject
{
    CentrosomeVis(std::shared_ptr<Centrosome> pCentrosome, GPUQueue* pQueue);

    virtual GPUMeshNode updateAndGetMeshNode() override;

private:
    void createCentrosomeGeometry();
    std::shared_ptr<Centrosome> m_pCentrosome;
    std::shared_ptr<GPUMesh> m_pUnitCylinderGpuMesh;
    // Cached root scene node; children[0] = X-axis node, children[1] = Y-axis node
    GPUMeshNode m_rootNode;

    // Build a linear transform that maps a unit cylinder (local Z = axis, local X/Y = radial)
    // to a cylinder aligned with 'axisUnit' in world space, with the given 'length' and 'radius'.
    // Returns a float3x3 whose rows are [radial1*radius, radial2*radius, axis*length].
    static float3x3 buildScaledCylinderMatrix(const float3& axisUnit, float length, float radius);

    // Update transforms for the two centriole cylinders (children[0], children[1])
    void updateCentriolesNodes();

    // Update/add child cylinders for each ring complex (children[2..])
    void updateRingComplexNodes();
};

