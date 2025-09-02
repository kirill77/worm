#pragma once

#include <string>
#include <memory>
#include "StringDict.h"

class TRNA
{
private:
    StringDict::ID m_id;         // StringDict ID that determines tRNA type (amino acid + anticodon)
    double m_fNumber;            // How much of this tRNA is available
    bool m_bCharged;             // Whether tRNA is loaded with amino acid
    double m_fChargingRate;      // Rate at which this tRNA gets charged with amino acid

public:
    TRNA(StringDict::ID id, double number, double chargingRate = 1.0)
        : m_id(id), m_fNumber(number), m_bCharged(false),
          m_fChargingRate(chargingRate) {}

    // Getters
    StringDict::ID getID() const { return m_id; }
    std::string getName() const { return StringDict::idToString(m_id); }
    std::string getAnticodon() const;
    double getNumber() const { return m_fNumber; }
    bool isCharged() const { return m_bCharged; }

    // tRNA charging
    void charge(double dt);
    void discharge();

    // Check if this tRNA matches a codon
    bool matchesCodon(const std::string& codon) const;
};

