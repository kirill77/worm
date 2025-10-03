#pragma once

#include <limits>
#include "vector.h"
#include "box.h"

// Slab-based ray–AABB intersection. Returns true if the infinite ray intersects the box.
// Outputs parametric distances along the ray (tNear, tFar) where pos + t * dir hits the box slabs.
// Constraints: dir components near zero are handled; if ray is parallel to an axis and outside
// the slab on that axis the function returns false.
template <class T, int n>
inline bool intersectRayAABB(const vector<T, n>& pos,
                             const vector<T, n>& dir,
                             const box<T, n>& b,
                             T& outMin,
                             T& outMax)
{
    const T negInf = std::numeric_limits<T>::lowest();
    const T posInf = std::numeric_limits<T>::max();
    T tNear = negInf;
    T tFar  = posInf;

    for (int i = 0; i < n; ++i)
    {
        const T o = pos[i];
        const T d = dir[i];
        const T minv = b.m_mins[i];
        const T maxv = b.m_maxs[i];

        if (std::abs(d) < std::numeric_limits<T>::epsilon())
        {
            // Parallel to slab; must be within bounds
            if (o < minv || o > maxv)
                return false;
            // No update to tNear/tFar
            continue;
        }

        const T invD = T(1) / d;
        T t1 = (minv - o) * invD;
        T t2 = (maxv - o) * invD;
        if (t1 > t2) std::swap(t1, t2);

        if (t1 > tNear) tNear = t1;
        if (t2 < tFar)  tFar  = t2;
        if (tNear > tFar) return false;
        if (tFar < T(0))  return false; // box is behind the ray
    }

    outMin = tNear;
    outMax = tFar;
    return true;
}


