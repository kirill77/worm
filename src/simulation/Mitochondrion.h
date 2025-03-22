#pragma once

#include "Organelle.h"
#include "Medium.h"
#include "math/vector.h"
#include <random>

class Mitochondrion : public Organelle
{
private:
    double m_fNumber;  // Number of mitochondria in the population
    
    static constexpr double ATP_PRODUCTION_RATE = 1.5;
    static constexpr double N_INITIAL_MITOCHONDRIA = 1e5;
    
    // Random number generator for ATP distribution
    std::mt19937 m_rng;

    // Generate ATP at n random positions in the medium
    void generateATP(Medium& medium, uint32_t n, double amount);
    // Helper method to generate a random position within the medium
    float3 generateRandomPosition();

public:
    Mitochondrion()
        : m_fNumber(N_INITIAL_MITOCHONDRIA) 
    {
        // Initialize random number generator with a seed
        std::random_device rd;
        m_rng.seed(rd());
    }

    void update(double dt, Cell& cell, Medium& medium) override;
};

