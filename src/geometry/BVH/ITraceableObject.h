#pragma once

#include <memory>
#include <cstdint>
#include "geometry/vectors/box.h"

// Forward declarations
struct ITraceableObject;
struct IRay;

struct IRay
{
    float3 m_vPos;
    float3 m_vDir;
    float m_fMin, m_fMax;
    virtual void notifyIntersection(float fDist, const ITraceableObject *pObject, uint32_t uSubObj) = 0;
};

// Can represent triangular mesh for example where every
// sub-object is a triangle. Must have at least one sub-object
struct ITraceableObject : std::enable_shared_from_this<ITraceableObject>
{
    uint32_t m_nSubObjects = 0;
    
    // bounding box of the whole object
    virtual box3 getBox() const = 0;
    
    // bounding box of a sub-object
    virtual box3 getSubObjectBox(uint32_t uSubObj) const = 0;
    
    virtual void trace(IRay& ray, uint32_t uSubObj) const = 0;
    
    virtual ~ITraceableObject() = default;
}; 