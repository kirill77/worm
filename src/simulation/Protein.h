#pragma once

#include <string>

struct Protein
{
    std::string m_sName;  // Name/type of the protein
    
    Protein(const std::string& name) : m_sName(name) {}
};

struct ProteinPopulation : public Protein
{
    double m_fNumber; // Number of molecules in this population
    
    ProteinPopulation(const std::string& sName, double fNumber)
        : Protein(sName), m_fNumber(fNumber) {}
};
