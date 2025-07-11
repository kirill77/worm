#pragma once

#include <memory>
#include "geometry/vectors/vector.h"

struct Centrosome;

// represents Gamma-TuRC (gamma-tubulin ring complex) - the place where a microtubule nucleates
struct Y_TuRC
{
    Y_TuRC(std::weak_ptr<Centrosome> pCentrosome);

    const float3& getDirection() { return m_vDir; }
    const float3& getPosition() { return m_vPosMicroM; }

private:
    std::weak_ptr<Centrosome> m_pCentrosome;

    float3 m_vDir; // normalized direction
    float3 m_vPosMicroM; // position in micro-meters in respect to centrosome center

    // microtubule grows by attaching alpha and beta tubulins
    uint32_t m_nAlphaTubulins = 0;
    uint32_t m_nBetaTubulins = 0;
};

