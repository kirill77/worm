#pragma once

#include <memory>
#include <vector>
#include <cassert>
#include "geometry/vectors/vector.h"
#include "geometry/mesh/MeshLocation.h"

// Physics representation of a microtubule as a sequence of points in 3D space
class PhysMicrotubule : public std::enable_shared_from_this<PhysMicrotubule>
{
public:
    // Microtubule dynamic state
    enum class MTState { Growing, Shrinking, Bound };

    PhysMicrotubule() = default;
    virtual ~PhysMicrotubule() = default;

    explicit PhysMicrotubule(std::vector<float3> points)
        : m_points(std::move(points))
    {
    }

    // Geometry accessors (const)
    const float3& getOrigin() const { return m_points[0]; }
    bool hasActiveMT() const { return m_points.size() >= 2; }
    float3 getTipPosition() const { return m_points.back(); }

    // Geometry accessors (non-const, for derived classes)
    std::vector<float3>& getPoints() { return m_points; }

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
    // Points defining the microtubule path (from minus end to plus end)
    std::vector<float3> m_points;
    
    // Current state (relevant for physics: Bound state enables cortical forces)
    MTState m_mtState = MTState::Growing;

    // Cortical attachment location (valid only when m_mtState == Bound)
    MeshLocation m_attachmentLocation;
};

