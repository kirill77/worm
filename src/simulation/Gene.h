#pragma once

#include <string>
#include <memory>
#include "MRNA.h"

class Gene
{
private:
    std::string m_sName;
    double m_fExpressionRate; // Rate of transcription
    double m_fBasalLevel;     // Basal expression level

public:
    Gene(const std::string& name, double expressionRate = 1.0, double basalLevel = 0.1)
        : m_sName(name),
          m_fExpressionRate(expressionRate), m_fBasalLevel(basalLevel) {}

    // Getters
    std::string getName() const { return m_sName; }
    double getExpressionRate() const { return m_fExpressionRate; }

    // Setters
    void setExpressionRate(double rate) { m_fExpressionRate = rate; }

    // Transcription
    std::shared_ptr<MRNA> transcribe(double dt) const;
};

