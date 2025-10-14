#pragma once

#include <memory>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <limits>
#include "Medium.h"
#include "CortexLocation.h"
#include "Organelle.h"
#include "geometry/vectors/vector.h"
#include "geometry/BVH/ITraceableObject.h"

/**
 * The Cortex class represents the cell membrane that separates
 * the internal cellular environment from the external environment.
 * It mediates interactions between internal and external media.
 */
class Cortex : public Organelle
{
public:


private:
    double m_fThickness; // Membrane thickness in micrometers
    std::shared_ptr<class BVHMesh> m_pCortexBVH;
    std::shared_ptr<class EdgeMesh> m_pCortexMesh;
    std::vector<CortexMolecules> m_pBindingSites;
    std::vector<Molecule> m_bindableMolecules;

public:
    /**
     * Constructor that initializes the cortex.
     *
     * @param pCell Weak pointer to the cell containing this cortex
     * @param fThickness Membrane thickness in micrometers
     */
    Cortex(std::weak_ptr<Cell> pCell, double fThickness = 0.01); // Default 10nm thickness

    /**
     * Update the cortex state.
     * This method updates cortex dynamics and can be extended to include
     * passive transport, signal transduction, etc.
     *
     * @param fDtSec Time step in seconds
     * @param cell Reference to the cell containing this cortex
     */
    void update(double fDtSec, Cell& cell) override;

    // Provide cortex surface mesh from external simulation (e.g., TensionSphere)
    void setMesh(const std::shared_ptr<EdgeMesh> pMesh);

    // Access underlying cortex surface mesh
    std::shared_ptr<class EdgeMesh> getEdgeMesh() const;

    /**
     * Initialize binding sites in the cell's internal medium.
     * This creates binding sites throughout the medium
     * that allow proteins to bind to the cell membrane surface.
     *
     * @param totalAmount Total amount of binding sites to distribute
     * @return True if binding sites were successfully added
     */
    bool initializeBindingSites(double totalAmount = 1000.0);

    /**
     * Transfer all molecules from cortex binding sites into the cell's internal medium grid
     * at the corresponding surface positions. After transfer, binding site populations are cleared.
     */
    void transferBindingSiteMoleculesToMedium();

    // Map normalized coordinates [-1,1] to cortex surface world position via ray cast
    float3 normalizedToWorld(const float3& normalizedPos);

    // Map world position to normalized coordinates [-1,1] using cortex bounding box
    float3 worldToNormalized(const float3& worldPos, bool isOnCortex = false) const;

    // Expose BVH mesh for visualization
    std::shared_ptr<class BVHMesh> getBVHMesh() const { return m_pCortexBVH; }

    // Ray used for cortex BVH tracing with intersection data
    struct CortexRay : public IRay {
        float distance = std::numeric_limits<float>::max();  // Distance to closest intersection
        uint32_t triangleIndex = 0;   // Index of intersected triangle
        float3 worldHitPoint = float3(0, 0, 0);  // World position of intersection
        bool hasHit = false;          // Whether intersection occurred

        // Constructor to set up ray parameters
        CortexRay(const float3& origin, const float3& direction);

        // IRay interface implementation
        void notifyIntersection(float fDist, const ITraceableObject*, uint32_t uSubObj) override;

        // Convenience method to get distance (0 if no hit for backward compatibility)
        float getDistance() const { return hasHit ? distance : 0.0f; }
    };

    // Find closest intersection with cortex surface along a ray
    bool findClosestIntersection(CortexRay& ray) const;

private:
    // Convert triangle index and barycentric coordinates to normalized [-1,1] coordinates
    float3 baryToNormalized(uint32_t triangleIndex, const float3& barycentric) const;

    // Pull molecules from grid cells at binding-site positions into binding sites
    void pullBindingSiteMoleculesFromMedium();
};