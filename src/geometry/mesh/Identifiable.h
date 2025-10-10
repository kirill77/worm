#pragma once

#include <atomic>
#include <cstdint>
#include <memory>

// Reusable base that provides a unique ID and shared_from_this via the base type.
// Note: shared_from_this() will return std::shared_ptr<Identifiable>.
class Identifiable : public std::enable_shared_from_this<Identifiable>
{
public:
    Identifiable()
        : m_id(++s_nextId)
    {
    }

    // Non-copyable to avoid duplicate IDs; allow move semantics
    Identifiable(const Identifiable&) = delete;
    Identifiable& operator=(const Identifiable&) = delete;
    Identifiable(Identifiable&&) = default;
    Identifiable& operator=(Identifiable&&) = default;

    uint64_t getId() const noexcept { return m_id; }

private:
    uint64_t m_id;
    static std::atomic<uint64_t> s_nextId;
};


