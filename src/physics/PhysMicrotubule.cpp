#include "PhysMicrotubule.h"
#include "chemistry/molecules/simConstants.h"
#include <cmath>

float PhysMicrotubule::getLastSegmentLength() const
{
    auto count = getVertexCount();
    float3 diff = getVertexPosition(count - 1) - getVertexPosition(count - 2);
    return sqrtf(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);
}

float PhysMicrotubule::getMTLengthMicroM() const
{
    auto count = getVertexCount();
    const float segmentLength = static_cast<float>(MoleculeConstants::MT_SEGMENT_LENGTH_MICROM);
    // All segments except the last are full length
    float totalLength = segmentLength * static_cast<float>(count - 2);
    // Add the last segment
    totalLength += getLastSegmentLength();
    return totalLength;
}

