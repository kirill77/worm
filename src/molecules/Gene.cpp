#include "pch.h"
#include "Gene.h"
#include "MRNA.h"
#include "GeneWiki.h"
#include <random>

std::shared_ptr<MRNA> Gene::transcribe(double dt) const
{
    // Calculate amount of mRNA produced based on expression rate and time step
    double mRNAAmount = m_fExpressionRate * dt + m_fBasalLevel;

    // Add some stochastic variation (noise in gene expression)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> noise(1.0, 0.1); // 10% noise
    mRNAAmount *= noise(gen);

    // Create new mRNA with typical half-life of 2.0 time units
    return std::make_shared<MRNA>(
        m_sName,
        mRNAAmount,
        2.0,  // half-life
        1.0   // translation rate
    );
} 