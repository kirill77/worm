#pragma once

#include <vector>
#include <memory>
#include "geometry/vectors/box.h"
#include "ITraceableObject.h"

// inheriting BVH from ITraceableObject enables hierarchy of BVHs
class BVH : public ITraceableObject
{
public:
    BVH();
    
    std::vector<std::shared_ptr<ITraceableObject>>& accessObjects();

    // Build/rebuild the BVH hierarchy
    void rebuildHierarchy();
    
    // ITraceableObject interface
    void trace(IRay& ray, uint32_t uSubObj) override;
    virtual box3 getBox() override;
    virtual box3 getSubObjectBox(uint32_t uSubObj) override;

private:
    struct SubObj
    {
        ITraceableObject* pObj;
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

    std::vector<std::shared_ptr<ITraceableObject>> m_pObjects;
    std::unique_ptr<Node> m_pRoot;

    // Helper methods
    bool rayIntersectsBox(const IRay& ray, const box3& box);
    std::unique_ptr<Node> buildNode(std::vector<SubObj>& subObjects, int depth = 0);
    void traceNode(IRay& ray, const Node* node);
    box3 calculateBoundingBox(const std::vector<SubObj>& subObjects);
    int getLongestAxis(const box3& box);
};

