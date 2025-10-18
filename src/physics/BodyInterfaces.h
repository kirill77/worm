#pragma once

#include <cstddef>
#include <utility>
#include <vector>
#include "geometry/vectors/vector.h"

// Body-agnostic node/edge/face views to support different topologies

struct INodeView
{
    virtual ~INodeView() = default;
    virtual size_t size() const = 0;
    virtual double3 getPosition(uint32_t i) const = 0;
    virtual double3 getVelocity(uint32_t i) const = 0;
    virtual double  getMass(uint32_t i) const = 0;
    virtual void    setPosition(uint32_t i, const double3& p) = 0;
    virtual void    setVelocity(uint32_t i, const double3& v) = 0;
    virtual void    addForce(uint32_t i, const double3& f) = 0;
};

struct IEdgeView
{
    virtual ~IEdgeView() = default;
    virtual size_t edgeCount() const = 0;
    virtual std::pair<uint32_t,uint32_t> edge(uint32_t e) const = 0;
    virtual double restLength(uint32_t e) const = 0;
};

struct IFaceView
{
    virtual ~IFaceView() = default;
    virtual size_t faceCount() const = 0;
    virtual uint3  face(uint32_t f) const = 0;
};

// Null-object views to avoid pointer/null checks
struct NullEdgeView final : public IEdgeView
{
    size_t edgeCount() const override { return 0; }
    std::pair<uint32_t,uint32_t> edge(uint32_t) const override { return {0u,0u}; }
    double restLength(uint32_t) const override { return 0.0; }
};

struct NullFaceView final : public IFaceView
{
    size_t faceCount() const override { return 0; }
    uint3  face(uint32_t) const override { return uint3{0u,0u,0u}; }
};

struct IBody
{
    virtual ~IBody() = default;
    virtual INodeView& nodes() = 0;
    virtual const INodeView& nodes() const = 0;
    // Default implementations return null-object views
    virtual const IEdgeView& edges() const { static const NullEdgeView kNull; return kNull; }
    virtual IEdgeView& edges() { static NullEdgeView kNull; return kNull; }
    virtual const IFaceView& faces() const { static const NullFaceView kNull; return kNull; }
    virtual IFaceView& faces() { static NullFaceView kNull; return kNull; }
};


