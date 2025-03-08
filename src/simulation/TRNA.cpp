#include "pch.h"
#include "TRNA.h"
#include <algorithm>

void TRNA::charge(double dt)
{
    if (!m_bCharged)
    {
        // Probability of charging in this time step
        double chargeProbability = m_fChargingRate * dt;
        if ((double)rand() / RAND_MAX < chargeProbability)
        {
            m_bCharged = true;
        }
    }
}

void TRNA::discharge()
{
    m_bCharged = false;
}

bool TRNA::matchesCodon(const std::string& codon) const
{
    // Convert codon to anticodon (simplified - in reality this would be more complex)
    std::string complementary = codon;
    for (char& c : complementary)
    {
        switch (c)
        {
            case 'A': c = 'U'; break;
            case 'U': c = 'A'; break;
            case 'G': c = 'C'; break;
            case 'C': c = 'G'; break;
        }
    }
    return complementary == m_sAnticodon;
}
