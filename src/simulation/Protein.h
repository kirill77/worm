#pragma once

#include <string>

struct Protein
{
    std::string m_sName;  // Name/type of the protein
};

struct ProteinPopulation : public Protein
{
    double m_fNumber; // Number of molecules in this population
};
