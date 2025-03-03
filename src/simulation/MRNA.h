#pragma once

#include <string>

class MRNA
{
    std::string m_sGeneName; // Gene that produced this mRNA
    std::string m_sProteinName; // Protein it produces
    double m_fNumber; // Amount of this mRNA in the cell
    double m_fHalfLife; // How quickly it degrades (in time units)

public:
    MRNA(std::string sGeneName, double fNumber, double fHalfLife)
        : m_sGeneName(sGeneName), m_fNumber(fNumber), m_fHalfLife(fHalfLife)
    {
    }

    void degrade(double dt)
    {
        // Simple exponential decay model for mRNA degradation
        m_fNumber *= exp(-dt / m_fHalfLife);
    }

    std::string getGeneName() const { return m_sGeneName; }
};

