#pragma once

#include <memory>
#include <vector>
#include <cassert>
#include "geometry/vectors/vector.h"
#include "geometry/mesh/MeshLocation.h"
#include "geometry/mesh/Vertices.h"

// Physics representation of a microtubule as a 1D mesh (polyline)
class PhysMicrotubule : public Vertices
{
public:
    // Microtubule dynamic state
    enum class MTState { Growing, Shrinking, Bound };

    PhysMicrotubule() = default;
    virtual ~PhysMicrotubule() = default;

    // Microtubule-specific geometry accessors
    float3 getOrigin() const { return getVertexPosition(0); }
    float3 getTipPosition() const { return getVertexPosition(getVertexCount() - 1); }

    // State accessors
    MTState getState() const { return m_mtState; }
    void setState(MTState state) { m_mtState = state; }

    // Geometry calculations
    float getLastSegmentLength() const;
    float getMTLengthMicroM() const;

    // Cortical attachment accessors (only valid when bound to cortex)
    const MeshLocation& getAttachmentLocation() const 
    { 
        assert(m_mtState == MTState::Bound && "Attachment location only valid when microtubule is bound");
        return m_attachmentLocation; 
    }
    MeshLocation& getAttachmentLocation() 
    { 
        assert(m_mtState == MTState::Bound && "Attachment location only valid when microtubule is bound");
        return m_attachmentLocation; 
    }
    void setAttachmentLocation(const MeshLocation& location) 
    { 
        assert(m_mtState == MTState::Bound && "Must set state to Bound before setting attachment location");
        m_attachmentLocation = location; 
    }

private:
    // Current state (relevant for physics: Bound state enables cortical forces)
    MTState m_mtState = MTState::Growing;

    // Cortical attachment location (valid only when m_mtState == Bound)
    MeshLocation m_attachmentLocation;
};

