#pragma once

#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <memory>
#include "../math/vector.h"
#include "connectedMesh/connectedMesh.h"

/**
 * @brief Class for simulating tension forces in a spherical cell cortex using a geodesic sphere model
 * implemented with Winged-Edge Data Structure for efficient mesh traversal
 */
class TensionSphere {
public:
    /**
     * @brief Constructor
     * @param subdivisionLevel Number of times to subdivide the base icosahedron (increases detail)
     */
    TensionSphere(uint32_t subdivisionLevel = 2);

    void makeTimeStep(double fDtSec);

    /**
     * @brief Get the underlying ConnectedMesh
     * @return Shared pointer to the ConnectedMesh
     */
    std::shared_ptr<ConnectedMesh> getConnectedMesh() const { return m_pMesh; }

private:
    // The underlying mesh data structure
    std::shared_ptr<ConnectedMesh> m_pMesh;
};

