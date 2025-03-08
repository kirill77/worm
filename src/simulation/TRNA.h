#pragma once

#include <string>
#include <memory>

class TRNA
{
private:
    std::string m_sAminoAcid;    // Amino acid carried by this tRNA
    std::string m_sAnticodon;    // Anticodon sequence
    double m_fNumber;            // How much of this tRNA is available
    bool m_bCharged;             // Whether tRNA is loaded with amino acid
    double m_fChargingRate;      // Rate at which this tRNA gets charged with amino acid

public:
    TRNA(const std::string& aminoAcid, const std::string& anticodon,
         double number, double chargingRate = 1.0)
        : m_sAminoAcid(aminoAcid), m_sAnticodon(anticodon),
          m_fNumber(number), m_bCharged(false),
          m_fChargingRate(chargingRate) {}

    // Getters
    std::string getAminoAcid() const { return m_sAminoAcid; }
    std::string getAnticodon() const { return m_sAnticodon; }
    double getNumber() const { return m_fNumber; }
    bool isCharged() const { return m_bCharged; }

    // tRNA charging
    void charge(double dt);
    void discharge();

    // Check if this tRNA matches a codon
    bool matchesCodon(const std::string& codon) const;
};

