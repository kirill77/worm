#include "BVH.h"
#include <algorithm>

std::vector<BVH::IObject*>& BVH::accessObjects()
{
    return m_pObjects;
}

void BVH::trace(IRay& ray)
{
    if (m_pRoot == nullptr)
    {
        // Fallback to linear traversal if no hierarchy is built
        for (IObject* pObject : m_pObjects)
        {
            if (pObject != nullptr)
            {
                box3 objectBox = pObject->getBox();
                if (rayIntersectsBox(ray, objectBox))
                {
                    pObject->trace(ray);
                }
            }
        }
    }
    else
    {
        // Use hierarchical traversal
        traceNode(ray, m_pRoot.get());
    }
}

void BVH::rebuildHierarchy()
{
    if (m_pObjects.empty())
    {
        m_pRoot = nullptr;
        return;
    }

    // Create a copy of the objects vector for building
    std::vector<IObject*> objects = m_pObjects;

    // Build the tree recursively
    m_pRoot = buildNode(objects, 0);
}

bool BVH::rayIntersectsBox(const IRay& ray, const box3& box)
{
    // Ray-box intersection test using slab method
    float3 invDir = float3(1.0f / ray.m_vDir.x, 1.0f / ray.m_vDir.y, 1.0f / ray.m_vDir.z);
    
    float3 t1 = (box.m_mins - ray.m_vPos) * invDir;
    float3 t2 = (box.m_maxs - ray.m_vPos) * invDir;
    
    float3 tMin = min(t1, t2);
    float3 tMax = max(t1, t2);
    
    // Manually compute max and min components to avoid template conflicts
    float tNear = std::max(tMin.x, std::max(tMin.y, tMin.z));
    float tFar = std::min(tMax.x, std::min(tMax.y, tMax.z));
    
    // Check if intersection is within ray's active range
    return tNear <= tFar && tFar >= ray.m_fMin && tNear <= ray.m_fMax;
}

std::unique_ptr<BVH::Node> BVH::buildNode(std::vector<IObject*>& objects, int depth)
{
    auto node = std::make_unique<Node>();
    
    // Calculate bounding box for this node
    node->m_boundingBox = calculateBoundingBox(objects);
    
    // Leaf node criteria: few objects or max depth reached
    const int MAX_LEAF_OBJECTS = 4;
    const int MAX_DEPTH = 20;
    
    if (objects.size() <= MAX_LEAF_OBJECTS || depth >= MAX_DEPTH)
    {
        // Create leaf node
        node->m_objects = objects;
        return node;
    }
    
    // Split objects along longest axis
    int axis = getLongestAxis(node->m_boundingBox);
    
    // Sort objects by their center position along the split axis
    std::sort(objects.begin(), objects.end(), [axis](IObject* a, IObject* b) {
        box3 boxA = a->getBox();
        box3 boxB = b->getBox();
        float3 centerA = (boxA.m_mins + boxA.m_maxs) * 0.5f;
        float3 centerB = (boxB.m_mins + boxB.m_maxs) * 0.5f;
        return centerA[axis] < centerB[axis];
    });
    
    // Split at median
    size_t mid = objects.size() / 2;
    std::vector<IObject*> leftObjects(objects.begin(), objects.begin() + mid);
    std::vector<IObject*> rightObjects(objects.begin() + mid, objects.end());
    
    // Recursively build child nodes
    node->m_pLeft = buildNode(leftObjects, depth + 1);
    node->m_pRight = buildNode(rightObjects, depth + 1);
    
    return node;
}

void BVH::traceNode(IRay& ray, const Node* node)
{
    if (node == nullptr)
        return;
        
    // Check if ray intersects this node's bounding box
    if (!rayIntersectsBox(ray, node->m_boundingBox))
        return;
    
    if (node->isLeaf())
    {
        // Leaf node - test all objects
        for (IObject* pObject : node->m_objects)
        {
            if (pObject != nullptr)
            {
                pObject->trace(ray);
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

box3 BVH::calculateBoundingBox(const std::vector<IObject*>& objects)
{
    if (objects.empty())
        return box3();
    
    box3 result = objects[0]->getBox();
    
    for (size_t i = 1; i < objects.size(); ++i)
    {
        if (objects[i] != nullptr)
        {
            result = result | objects[i]->getBox();
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
