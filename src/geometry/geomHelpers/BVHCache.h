#pragma once

#include <memory>
#include <unordered_map>
#include <mutex>
#include <cstdint>

// Simple flyweight-style cache for BVHMesh per TriangleMesh instance.
// Rebuilds BVHMesh whenever TriangleMesh version changes.
class TriangleMesh;
class BVHMesh;

class BVHCache
{
public:
    static BVHCache& instance();

    std::shared_ptr<BVHMesh> getOrCreate(const std::shared_ptr<TriangleMesh>& mesh);

private:
    BVHCache() = default;

    struct Entry
    {
        std::weak_ptr<TriangleMesh>    mesh;
        std::weak_ptr<BVHMesh> bvh;
        uint64_t               cachedVersion = 0;
    };

    std::unordered_map<uint64_t, Entry> m_entries;
    std::mutex m_mutex;
};


