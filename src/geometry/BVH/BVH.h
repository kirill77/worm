#pragma once

#include <vector>
#include <memory>
#include "geometry/vectors/box.h"

class BVH
{
    struct IObject;

    struct IRay
    {
        float3 m_vPos;
        float3 m_vDir;
        float m_fMin, m_fMax;
        void* m_pContext;
        virtual void notifyIntersection(float fDist, IObject* pObject) = 0;
    };
    struct IObject
    {
        virtual box3 getBox() = 0;
        virtual void trace(IRay& ray) = 0;
    };

    std::vector<IObject*>& accessObjects();

    void trace(IRay& ray);

    // Build/rebuild the BVH hierarchy
    void rebuildHierarchy();

private:
    struct Node
    {
        box3 m_boundingBox;
        std::unique_ptr<Node> m_pLeft;
        std::unique_ptr<Node> m_pRight;
        std::vector<IObject*> m_objects;  // Only used for leaf nodes
        
        bool isLeaf() const { return m_pLeft == nullptr && m_pRight == nullptr; }
    };

    std::vector<IObject*> m_pObjects;
    std::unique_ptr<Node> m_pRoot;

    // Helper methods
    bool rayIntersectsBox(const IRay& ray, const box3& box);
    std::unique_ptr<Node> buildNode(std::vector<IObject*>& objects, int depth = 0);
    void traceNode(IRay& ray, const Node* node);
    box3 calculateBoundingBox(const std::vector<IObject*>& objects);
    int getLongestAxis(const box3& box);
};

