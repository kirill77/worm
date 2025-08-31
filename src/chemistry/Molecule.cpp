#include "Molecule.h"
#include "StringDict.h"

// Constructor with name and type - automatically tries to use ID optimization
Molecule::Molecule(const std::string& name, ChemicalType type) : m_type(type)
{
    // Try to convert string to StringDict ID first for optimization
    m_id = StringDict::stringToId(name);

    assert(m_type != ChemicalType::OTHER);

    if (m_id != StringDict::ID::eUNKNOWN) {
        // Found matching ID - use optimized storage (empty string)
        m_sName.clear();  // Explicitly clear for safety
    } else {
        // Unknown molecule - store the string
        m_sName = name;
    }
}