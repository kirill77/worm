#include "BVH.h"
#include <algorithm>
#include "geometry/vectors/intersections.h"

BVH::BVH()
{
    // BVH has exactly one sub-object (itself)
    m_nSubObjects = 1;
}

std::vector<std::shared_ptr<ITraceableObject>>& BVH::accessObjects()
{
    return m_pObjects;
}

box3 BVH::getBox() const
{
    assert(m_pRoot != nullptr);
    return m_pRoot->m_boundingBox;
}

box3 BVH::getSubObjectBox(uint32_t uSubObj) const
{
    // BVH has only one sub-object (itself)
    assert(uSubObj == 0);
    return getBox();
}

void BVH::trace(IRay& ray, uint32_t uSubObj) const
{
    assert(m_pRoot != nullptr && uSubObj == 0);
    // Use hierarchical traversal
    traceNode(ray, m_pRoot.get());
}

void BVH::rebuildHierarchy()
{
    if (m_pObjects.empty())
    {
        m_pRoot = nullptr;
        return;
    }

    // Expand all objects into their sub-objects
    std::vector<SubObj> subObjects;
    
    // Calculate total number of sub-objects to reserve space
    size_t totalSubObjects = 0;
    for (auto pObject : m_pObjects)
    {
        if (pObject != nullptr && pObject->m_nSubObjects > 0)
        {
            totalSubObjects += pObject->m_nSubObjects;
        }
    }
    subObjects.reserve(totalSubObjects);
    
    for (auto pObject : m_pObjects)
    {
        if (pObject != nullptr && pObject->m_nSubObjects > 0)
        {
            for (uint32_t i = 0; i < pObject->m_nSubObjects; ++i)
            {
                SubObj subObj;
                subObj.pObj = pObject.get();
                subObj.m_uSubObj = i;
                subObjects.push_back(subObj);
            }
        }
    }

    if (subObjects.empty())
    {
        m_pRoot = nullptr;
        return;
    }

    // Build the tree recursively
    m_pRoot = buildNode(subObjects, 0);
}

bool BVH::rayIntersectsBox(const IRay& ray, const box3& box) const
{
    float tNear, tFar;
    if (!intersectRayAABB<float, 3>(ray.m_vPos, ray.m_vDir, box, tNear, tFar))
        return false;
    // Check if intersection interval overlaps the ray's active range
    return tFar >= ray.m_fMin && tNear <= ray.m_fMax;
}

std::unique_ptr<BVH::Node> BVH::buildNode(std::vector<SubObj>& subObjects, int depth)
{
    auto node = std::make_unique<Node>();
    
    // Calculate bounding box for this node
    node->m_boundingBox = calculateBoundingBox(subObjects);
    
    // Leaf node criteria: few objects or max depth reached
    const int MAX_LEAF_OBJECTS = 4;
    const int MAX_DEPTH = 20;
    
    if (subObjects.size() <= MAX_LEAF_OBJECTS || depth >= MAX_DEPTH)
    {
        // Create leaf node
        node->m_pSubObjects = subObjects;
        return node;
    }
    
    // Split objects along longest axis
    int axis = getLongestAxis(node->m_boundingBox);
    
    // Sort sub-objects by their center position along the split axis
    std::sort(subObjects.begin(), subObjects.end(), [axis](const SubObj& a, const SubObj& b) {
        box3 boxA = a.pObj->getSubObjectBox(a.m_uSubObj);
        box3 boxB = b.pObj->getSubObjectBox(b.m_uSubObj);
        float3 centerA = (boxA.m_mins + boxA.m_maxs) * 0.5f;
        float3 centerB = (boxB.m_mins + boxB.m_maxs) * 0.5f;
        return centerA[axis] < centerB[axis];
    });
    
    // Split at median
    size_t mid = subObjects.size() / 2;
    std::vector<SubObj> leftObjects(subObjects.begin(), subObjects.begin() + mid);
    std::vector<SubObj> rightObjects(subObjects.begin() + mid, subObjects.end());
    
    // Recursively build child nodes
    node->m_pLeft = buildNode(leftObjects, depth + 1);
    node->m_pRight = buildNode(rightObjects, depth + 1);
    
    return node;
}

void BVH::traceNode(IRay& ray, const Node* node) const
{
    if (node == nullptr)
        return;
        
    // Check if ray intersects this node's bounding box
    if (!rayIntersectsBox(ray, node->m_boundingBox))
        return;
    
    if (node->isLeaf())
    {
        // Leaf node - test all sub-objects
        for (const SubObj& subObj : node->m_pSubObjects)
        {
            if (subObj.pObj != nullptr)
            {
                // Test the bounding box of this specific sub-object first
                box3 subObjBox = subObj.pObj->getSubObjectBox(subObj.m_uSubObj);
                if (rayIntersectsBox(ray, subObjBox))
                {
                    // Let the object handle ray tracing for this specific sub-object
                    // The object's trace method should internally check which sub-objects to test
                    subObj.pObj->trace(ray, subObj.m_uSubObj);
                }
            }
        }
    }
    else
    {
        // Internal node - recursively traverse children
        traceNode(ray, node->m_pLeft.get());
        traceNode(ray, node->m_pRight.get());
    }
}

box3 BVH::calculateBoundingBox(const std::vector<SubObj>& subObjects)
{
    if (subObjects.empty())
        return box3();
    
    box3 result = subObjects[0].pObj->getSubObjectBox(subObjects[0].m_uSubObj);
    
    for (size_t i = 1; i < subObjects.size(); ++i)
    {
        if (subObjects[i].pObj != nullptr)
        {
            result = result | subObjects[i].pObj->getSubObjectBox(subObjects[i].m_uSubObj);
        }
    }
    
    return result;
}



int BVH::getLongestAxis(const box3& box)
{
    float3 extent = box.m_maxs - box.m_mins;
    
    if (extent.x >= extent.y && extent.x >= extent.z)
        return 0; // X axis
    else if (extent.y >= extent.z)
        return 1; // Y axis
    else
        return 2; // Z axis
}
