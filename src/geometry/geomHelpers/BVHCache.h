#pragma once

#include <memory>
#include <unordered_map>
#include <mutex>
#include <cstdint>

// Simple flyweight-style cache for BVHMesh per Mesh instance.
// Rebuilds BVHMesh whenever Mesh version changes.
class Mesh;
class BVHMesh;

class BVHCache
{
public:
    static BVHCache& instance();

    std::shared_ptr<BVHMesh> getOrCreate(const std::shared_ptr<Mesh>& mesh);

private:
    BVHCache() = default;

    struct Entry
    {
        std::weak_ptr<Mesh>    mesh;
        std::weak_ptr<BVHMesh> bvh;
        uint64_t               cachedVersion = 0;
    };

    std::unordered_map<uint64_t, Entry> m_entries;
    std::mutex m_mutex;
};


