#include "BVHCache.h"
#include "BVHMesh.h"
#include "geometry/mesh/Mesh.h"

BVHCache& BVHCache::instance()
{
    static BVHCache cache;
    return cache;
}

std::shared_ptr<BVHMesh> BVHCache::getOrCreate(const std::shared_ptr<Mesh>& mesh)
{
    if (!mesh) return nullptr;
    std::lock_guard<std::mutex> lock(m_mutex);

    const uint64_t key = mesh->getId();
    auto it = m_entries.find(key);
    const uint64_t ver = mesh->getVersion();

    if (it == m_entries.end())
    {
        auto bvh = std::make_shared<BVHMesh>(mesh);
        bvh->rebuildForCurrentMesh();
        m_entries.emplace(key, Entry{ std::weak_ptr<Mesh>(mesh), std::weak_ptr<BVHMesh>(bvh), ver });
        return bvh;
    }

    Entry& e = it->second;

    // If mesh expired or replaced, rebuild and replace entry
    if (e.mesh.expired())
    {
        auto bvh = std::make_shared<BVHMesh>(mesh);
        bvh->rebuildForCurrentMesh();
        e = Entry{ std::weak_ptr<Mesh>(mesh), std::weak_ptr<BVHMesh>(bvh), ver };
        return bvh;
    }

    auto bvhLocked = e.bvh.lock();
    if (!bvhLocked)
    {
        auto bvh = std::make_shared<BVHMesh>(mesh);
        bvh->rebuildForCurrentMesh();
        e.bvh = bvh;
        e.cachedVersion = ver;
        return bvh;
    }

    if (e.cachedVersion != ver)
    {
        // Mesh changed: rebuild BVHMesh
        auto bvh = std::make_shared<BVHMesh>(mesh);
        bvh->rebuildForCurrentMesh();
        e.bvh = bvh;
        e.cachedVersion = ver;
        return bvh;
    }

    return bvhLocked;
}


