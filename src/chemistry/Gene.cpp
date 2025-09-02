#include "pch.h"
#include "Gene.h"
#include "MoleculeWiki.h"
#include "GeneWiki.h"
#include "StringDict.h"
#include <random>

std::shared_ptr<MPopulation> Gene::transcribe(double dt) const
{
    // Calculate amount of mRNA produced based on expression rate and time step
    double mRNAAmount = m_fExpressionRate * dt + m_fBasalLevel;

    // Add some stochastic variation (noise in gene expression)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> noise(1.0, 0.1); // 10% noise
    mRNAAmount *= noise(gen);

    // Create new RNA molecule population
    Molecule rnaMolecule(m_id, ChemicalType::MRNA);
    return std::make_shared<MPopulation>(rnaMolecule, mRNAAmount);
} 