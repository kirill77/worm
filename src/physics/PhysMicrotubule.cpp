#include "PhysMicrotubule.h"
#include "chemistry/molecules/simConstants.h"
#include <cmath>

float PhysMicrotubule::getLastSegmentLength() const
{
    if (m_points.size() < 2) return 0.0f;
    float3 diff = m_points.back() - m_points[m_points.size() - 2];
    return sqrtf(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);
}

float PhysMicrotubule::getMTLengthMicroM() const
{
    if (m_points.size() < 2) return 0.0f;
    
    const float segmentLength = static_cast<float>(MoleculeConstants::MT_SEGMENT_LENGTH_MICROM);
    // All segments except the last are full length
    float totalLength = segmentLength * static_cast<float>(m_points.size() - 2);
    // Add the last segment
    totalLength += getLastSegmentLength();
    return totalLength;
}

