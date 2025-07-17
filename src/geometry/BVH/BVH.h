#pragma once

#include <vector>
#include <memory>
#include "geometry/vectors/box.h"

class BVH
{
public:
    struct IObject;

    struct IRay
    {
        float3 m_vPos;
        float3 m_vDir;
        float m_fMin, m_fMax;
        virtual void notifyIntersection(float fDist, IObject *pObject, uint32_t uSubObj) = 0;
    };
    
    // Can represent triangular mesh for example where every
    // sub-object is a triangle. Must have at least one sub-object
    struct IObject : std::enable_shared_from_this<IObject>
    {
        uint32_t m_nSubObjects = 0;
        // bounding box of the whole object
        virtual box3 getBox() = 0;
        // bounding box of a sub-object
        virtual box3 getSubObjectBox(uint32_t uSubObj) = 0;
        virtual void trace(IRay& ray) = 0;
    };

    std::vector<std::shared_ptr<IObject>>& accessObjects();

    void trace(IRay& ray);

    // Build/rebuild the BVH hierarchy
    void rebuildHierarchy();

private:
    struct SubObj
    {
        IObject* pObj;
        uint32_t m_uSubObj;
    };
    struct Node
    {
        box3 m_boundingBox;
        std::unique_ptr<Node> m_pLeft;
        std::unique_ptr<Node> m_pRight;
        std::vector<SubObj> m_pSubObjects;  // Only used for leaf nodes
        
        bool isLeaf() const { return m_pLeft == nullptr && m_pRight == nullptr; }
    };

    std::vector<std::shared_ptr<IObject>> m_pObjects;
    std::unique_ptr<Node> m_pRoot;

    // Helper methods
    bool rayIntersectsBox(const IRay& ray, const box3& box);
    std::unique_ptr<Node> buildNode(std::vector<SubObj>& subObjects, int depth = 0);
    void traceNode(IRay& ray, const Node* node);
    box3 calculateBoundingBox(const std::vector<SubObj>& subObjects);
    int getLongestAxis(const box3& box);
};

