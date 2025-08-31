#include "pch.h"
#include "Spindle.h"
#include "Cell.h"
#include "Medium.h"
#include "chemistry/StringDict.h"
#include <algorithm>

void Spindle::update(double dt, Cell& cell)
{
    Medium& medium = cell.getInternalMedium();
    
    // Update assembly state first
    updateAssembly(dt, cell, medium);
    
    // Only update position and orientation if assembled
    if (m_isAssembled)
    {
        updatePosition(dt, cell, medium);
        updateOrientation(dt, cell, medium);
    }
}

void Spindle::updateAssembly(double dt, Cell& cell, Medium& medium)
{
    // Only grow if not fully assembled and we have ATP
    if (!m_isAssembled && m_length < MAX_LENGTH)
    {
        float growthAmount = GROWTH_RATE * static_cast<float>(dt);
        float atpNeeded = growthAmount * ATP_PER_GROWTH;
        
        if (cell.consumeATP(atpNeeded))
        {
            m_length = std::min(MAX_LENGTH, m_length + growthAmount);
            
            // Check if assembly is complete
            if (m_length >= MAX_LENGTH)
            {
                m_isAssembled = true;
            }
        }
    }
}

void Spindle::updatePosition(double dt, Cell& cell, Medium& medium)
{
    // Calculate net force from cortical pulling
    float3 corticalForce = calculateCorticalForces(medium);
    
    // Update position based on force
    // Simple damped motion: force = -velocity (assuming mass = 1)
    float3 velocity = corticalForce * static_cast<float>(dt);
    m_position = m_position + velocity;
    
    // Keep within cell bounds (-1 to 1)
    for (int i = 0; i < 3; ++i)
    {
        m_position[i] = std::max(-0.9f, std::min(0.9f, m_position[i]));
    }
}

void Spindle::updateOrientation(double dt, Cell& cell, Medium& medium)
{
    // Get preferred orientation based on PAR polarity
    float3 preferredDir = calculatePreferredOrientation(medium);
    
    // Calculate rotation needed
    float3 rotationAxis = cross(m_orientation, preferredDir);
    float rotationAmount = length(rotationAxis);
    
    if (rotationAmount > 0.001f)  // Only rotate if significant difference
    {
        // Normalize rotation axis
        rotationAxis = rotationAxis * (1.0f / rotationAmount);
        
        // Calculate rotation step (with ATP cost)
        float rotationStep = ROTATION_RATE * static_cast<float>(dt);
        float atpNeeded = rotationStep * ATP_PER_ROTATION;
        
        if (cell.consumeATP(atpNeeded))
        {
            // Apply rotation
            float angle = std::min(rotationAmount, rotationStep);
            float cosAngle = cos(angle);
            float sinAngle = sin(angle);
            
            // Rodrigues rotation formula
            m_orientation = m_orientation * cosAngle + 
                          cross(rotationAxis, m_orientation) * sinAngle +
                          rotationAxis * dot(rotationAxis, m_orientation) * (1.0f - cosAngle);
            
            // Ensure normalized
            m_orientation = normalize(m_orientation);
        }
    }
}

float3 Spindle::calculateCorticalForces(const Medium& medium) const
{
    switch (m_cellType)
    {
        case CellType::Zygote:
        case CellType::Germline1:
        case CellType::Germline2:
        case CellType::Germline3:
            return calculatePARBasedForces(medium);
            
        default:
            return calculateDefaultForces(medium);
    }
}

float3 Spindle::calculatePARBasedForces(const Medium& medium) const
{
    float3 netForce(0.0f);
    
    // Get positions of spindle poles
    float3 minusPole = getMinusPole();
    float3 plusPole = getPlusPole();
    
    // Sample cortical points
    static const int NUM_SAMPLES = 8;
    for (int i = 0; i < NUM_SAMPLES; ++i)
    {
        float angle = 2.0f * 3.14159f * i / NUM_SAMPLES;
        
        // Sample anterior cortex
        float3 anteriorPoint(cos(angle), 0.95f, sin(angle));
        double anteriorPARs = medium.getMoleculeNumber(Molecule(StringDict::idToString(StringDict::ID::PAR_3), ChemicalType::PROTEIN), anteriorPoint) +
                              medium.getMoleculeNumber(Molecule(StringDict::idToString(StringDict::ID::PAR_6), ChemicalType::PROTEIN), anteriorPoint) +
                              medium.getMoleculeNumber(Molecule(StringDict::idToString(StringDict::ID::PKC_3), ChemicalType::PROTEIN), anteriorPoint);

        // Force on minus pole from anterior
        float3 toMinusPole = minusPole - anteriorPoint;
        float minusDist = length(toMinusPole);
        netForce = netForce + normalize(toMinusPole) * 
                  (float)(CORTICAL_FORCE * anteriorPARs / (minusDist * minusDist));
        
        // Sample posterior cortex
        float3 posteriorPoint(cos(angle), -0.95f, sin(angle));
        double posteriorPARs = medium.getMoleculeNumber(Molecule(StringDict::idToString(StringDict::ID::PAR_1), ChemicalType::PROTEIN), posteriorPoint) +
                               medium.getMoleculeNumber(Molecule(StringDict::idToString(StringDict::ID::PAR_2), ChemicalType::PROTEIN), posteriorPoint);

        // Force on plus pole from posterior (stronger in P lineage)
        float3 toPlusPole = plusPole - posteriorPoint;
        float plusDist = length(toPlusPole);
        double forceMagnitude = CORTICAL_FORCE * posteriorPARs;
        if (m_cellType == CellType::Zygote) forceMagnitude *= 1.5; // Stronger in zygote
        
        netForce = netForce + normalize(toPlusPole) * 
                  (float)(forceMagnitude / (plusDist * plusDist));
    }
    
    return netForce;
}

float3 Spindle::calculateDefaultForces(const Medium& medium) const
{
    float3 netForce(0.0f);
    
    // For somatic cells, implement symmetric pulling forces
    // This tends to center and orient the spindle along the long axis
    
    // Get positions of spindle poles
    float3 minusPole = getMinusPole();
    float3 plusPole = getPlusPole();
    
    // Sample cortical points
    static const int NUM_SAMPLES = 8;
    for (int i = 0; i < NUM_SAMPLES; ++i)
    {
        float angle = 2.0f * 3.14159f * i / NUM_SAMPLES;
        
        // Sample points around the cortex
        for (float y = -0.95f; y <= 0.95f; y += 1.9f)
        {
            float3 cortexPoint(cos(angle), y, sin(angle));
            
            // Equal forces on both poles
            float3 toMinusPole = minusPole - cortexPoint;
            float3 toPlusPole = plusPole - cortexPoint;
            
            float minusDist = length(toMinusPole);
            float plusDist = length(toPlusPole);
            
            netForce = netForce + normalize(toMinusPole) * (CORTICAL_FORCE / (minusDist * minusDist)) +
                                 normalize(toPlusPole) * (CORTICAL_FORCE / (plusDist * plusDist));
        }
    }
    
    return netForce;
}

float3 Spindle::calculatePreferredOrientation(const Medium& medium) const
{
    switch (m_cellType)
    {
        case CellType::Zygote:
        case CellType::Germline1:
        case CellType::Germline2:
        case CellType::Germline3:
        {
            // PAR-based orientation for germline lineage
            float3 anterior(0.0f, 0.95f, 0.0f);
            float3 posterior(0.0f, -0.95f, 0.0f);
            
            double anteriorPARs = medium.getMoleculeNumber(Molecule(StringDict::idToString(StringDict::ID::PAR_3), ChemicalType::PROTEIN), anterior) +
                                  medium.getMoleculeNumber(Molecule(StringDict::idToString(StringDict::ID::PAR_6), ChemicalType::PROTEIN), anterior) +
                                  medium.getMoleculeNumber(Molecule(StringDict::idToString(StringDict::ID::PKC_3), ChemicalType::PROTEIN), anterior);
            
            double posteriorPARs = medium.getMoleculeNumber(Molecule(StringDict::idToString(StringDict::ID::PAR_1), ChemicalType::PROTEIN), posterior) +
                                   medium.getMoleculeNumber(Molecule(StringDict::idToString(StringDict::ID::PAR_2), ChemicalType::PROTEIN), posterior);
            
            float3 direction = normalize(posterior - anterior);
            return direction * (posteriorPARs > anteriorPARs ? 1.0f : -1.0f);
        }
        
        default:
            // For somatic cells, orient along the cell's long axis
            return normalize(float3(0.0f, 1.0f, 0.0f));
    }
} 