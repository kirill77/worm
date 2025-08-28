#pragma once

#include <string>
#include <memory>
#include <assert.h>

// Forward declaration
class BindingSurface;

struct MPopulation
{
    std::string m_sName;  // Name/type of the protein
    double m_fNumber;     // Number of molecules in this population
    
    MPopulation(const std::string& sName, double fNumber)
        : m_sName(sName), m_fNumber(fNumber) {}
    
    // Check if this protein population is bound to any surface
    bool isBound() const
    {
        return !m_pBindingSurface.expired();
    }

    // Get a shared_ptr to the binding surface (if bound)
    std::shared_ptr<BindingSurface> getBindingSurface() const {
        return m_pBindingSurface.lock();
    }

    // Bind this protein population to a surface
    void bindTo(std::shared_ptr<BindingSurface> pSurface)
    {
        // if you want to bind to a different surface - have to unbind first
        assert(m_pBindingSurface.expired() || m_pBindingSurface.lock() == pSurface);
        m_pBindingSurface = pSurface;
    }

    // Unbind this protein population from its current surface
    void unbind()
    {
        m_pBindingSurface.reset();
    }

private:
    // Weak pointer to the surface this protein population is bound to (if any)
    std::weak_ptr<BindingSurface> m_pBindingSurface;
};
