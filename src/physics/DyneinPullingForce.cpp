#include "DyneinPullingForce.h"
#include "PhysCentrosome.h"
#include "PhysMicrotubule.h"
#include "PhysicsMesh.h"
#include "geometry/mesh/TriangleMesh.h"
#include <cmath>

DyneinPullingForce::DyneinPullingForce(
    PhysicsMesh& cortexBody,
    const std::vector<std::shared_ptr<PhysCentrosome>>& centrosomes,
    double pullingForcePerMT)
    : m_cortexBody(cortexBody)
    , m_centrosomes(centrosomes)
    , m_pullingForcePerMT(pullingForcePerMT)
{
}

void DyneinPullingForce::apply(double dt)
{
    // Iterate through all centrosomes and their microtubules
    for (const auto& pCentrosome : m_centrosomes)
    {
        if (!pCentrosome) continue;

        // Get centrosome position in cell space (transform from normalized to cell coordinates)
        const affine3& toNormalized = pCentrosome->getToNormalizedCell();
        const float3& centrosomePos = toNormalized.m_translation;

        // Iterate through microtubules
        for (const auto& pMT : pCentrosome->getMicrotubules())
        {
            if (!pMT) continue;

            // Only apply force if microtubule is bound to cortex
            if (pMT->getState() != PhysMicrotubule::MTState::Bound) continue;

            // Get microtubule tip position (in centrosome-local space)
            if (!pMT->hasActiveMT()) continue;
            
            float3 tipLocal = pMT->getTipPosition();
            
            // Transform tip position from centrosome-local to cell space
            // For now, assume simple translation (no rotation)
            float3 tipCell = float3(
                centrosomePos.x + tipLocal.x,
                centrosomePos.y + tipLocal.y,
                centrosomePos.z + tipLocal.z
            );

            // Calculate pulling force direction: from tip toward centrosome
            float3 forceDir = float3(
                centrosomePos.x - tipCell.x,
                centrosomePos.y - tipCell.y,
                centrosomePos.z - tipCell.z
            );
            
            float length = sqrtf(forceDir.x * forceDir.x + forceDir.y * forceDir.y + forceDir.z * forceDir.z);
            if (length > 1e-6f)
            {
                // Normalize direction
                forceDir = float3(forceDir.x / length, forceDir.y / length, forceDir.z / length);
                
                // Scale by force magnitude
                float3 force = float3(
                    forceDir.x * static_cast<float>(m_pullingForcePerMT),
                    forceDir.y * static_cast<float>(m_pullingForcePerMT),
                    forceDir.z * static_cast<float>(m_pullingForcePerMT)
                );

                       // Get attachment location (triangle and barycentric coordinates)
                       const auto& mesh = m_cortexBody.m_pMesh;
                       const MeshLocation& attachmentLoc = pMT->getAttachmentLocation();
                       
                       // Get triangle vertices
                       uint3 triangle = mesh->getTriangleVertices(attachmentLoc.m_triangleIndex);
                
                       // Distribute force to triangle vertices using barycentric weights
                       const float3& baryCoords = attachmentLoc.getBarycentric();
                       auto& v0 = m_cortexBody.getVertex(triangle.x);
                       auto& v1 = m_cortexBody.getVertex(triangle.y);
                       auto& v2 = m_cortexBody.getVertex(triangle.z);
                       
                       v0.m_vForce.x += force.x * baryCoords.x;
                       v0.m_vForce.y += force.y * baryCoords.x;
                       v0.m_vForce.z += force.z * baryCoords.x;
                       
                       v1.m_vForce.x += force.x * baryCoords.y;
                       v1.m_vForce.y += force.y * baryCoords.y;
                       v1.m_vForce.z += force.z * baryCoords.y;
                       
                       v2.m_vForce.x += force.x * baryCoords.z;
                       v2.m_vForce.y += force.y * baryCoords.z;
                       v2.m_vForce.z += force.z * baryCoords.z;
            }
        }
    }
}

