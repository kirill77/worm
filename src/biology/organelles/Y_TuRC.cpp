#include "Y_TuRC.h"
#include "Centrosome.h"
#include <random>
#include <cmath>
#include <numbers>

Y_TuRC::Y_TuRC(std::weak_ptr<Centrosome> pCentrosome)
    : m_pCentrosome(pCentrosome)
    , m_nAlphaTubulins(0)
    , m_nBetaTubulins(0)
{
    // Initialize random direction for microtubule nucleation
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dis(-1.0f, 1.0f);
    
    // Generate random normalized direction
    m_vDir = float3(dis(gen), dis(gen), dis(gen));
    float length = sqrtf(m_vDir.x * m_vDir.x + m_vDir.y * m_vDir.y + m_vDir.z * m_vDir.z);
    if (length > 0.0f) {
        m_vDir = float3(m_vDir.x / length, m_vDir.y / length, m_vDir.z / length);
    } else {
        m_vDir = float3(1.0f, 0.0f, 0.0f); // Default direction
    }
    
    // Position randomly within PCM radius from the centrosome using spherical distribution
    float pcmRadius = 0.5f;  // Default radius
    if (auto centrosome = pCentrosome.lock()) {
        pcmRadius = centrosome->getPCMRadius();
    }
    
    // Generate uniform point within sphere
    std::uniform_real_distribution<float> unitDis(0.0f, 1.0f);
    std::uniform_real_distribution<float> thetaDis(0.0f, 2.0f * std::numbers::pi_v<float>);
    std::uniform_real_distribution<float> phiDis(-1.0f, 1.0f);
    
    float r = pcmRadius * std::cbrt(unitDis(gen));  // Cube root for uniform volume distribution
    float theta = thetaDis(gen);  // Azimuthal angle
    float phi = std::acos(phiDis(gen));  // Polar angle
    
    // Convert spherical to Cartesian coordinates
    m_vPosMicroM = float3(
        r * std::sin(phi) * std::cos(theta),
        r * std::sin(phi) * std::sin(theta), 
        r * std::cos(phi)
    );
}
