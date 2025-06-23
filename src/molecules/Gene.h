#pragma once

#include <string>
#include <memory>
#include "MRNA.h"
#include "StringDict.h"

class Gene
{
private:
    StringDict::ID m_id;
    double m_fExpressionRate; // Rate of transcription
    double m_fBasalLevel;     // Basal expression level

public:
    Gene(StringDict::ID id, double expressionRate = 1.0, double basalLevel = 0.1)
        : m_id(id),
          m_fExpressionRate(expressionRate), m_fBasalLevel(basalLevel) {}

    // Getters
    std::string getName() const { return StringDict::idToString(m_id); }
    StringDict::ID getId() const { return m_id; }
    double getExpressionRate() const { return m_fExpressionRate; }

    // Setters
    void setExpressionRate(double rate) { m_fExpressionRate = rate; }

    // Transcription
    std::shared_ptr<MRNA> transcribe(double dt) const;
};

