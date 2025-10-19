#pragma once

#include <memory>
#include <vector>
#include "geometry/vectors/affine.h"

class PhysMicrotubule;

// Physics representation of a centrosome - microtubule organizing center
// This class contains no biological details, only geometric and organizational data
class PhysCentrosome
{
public:
    PhysCentrosome() = default;

    // Transform accessors
    const affine3& getToNormalizedCell() const { return m_toNormalizedCell; }
    affine3& getToNormalizedCell() { return m_toNormalizedCell; }
    void setToNormalizedCell(const affine3& transform) { m_toNormalizedCell = transform; }

    // Microtubule accessors
    const std::vector<std::shared_ptr<PhysMicrotubule>>& getMicrotubules() const { return m_pMicrotubules; }
    std::vector<std::shared_ptr<PhysMicrotubule>>& getMicrotubules() { return m_pMicrotubules; }

private:
    // Transform from centrosome-local space to normalized cell space (-1,1)
    affine3 m_toNormalizedCell;

    // Microtubules nucleated from this centrosome
    std::vector<std::shared_ptr<PhysMicrotubule>> m_pMicrotubules;
};

