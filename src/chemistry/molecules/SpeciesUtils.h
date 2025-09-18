#pragma once

#include <string>
#include <assert.h>
#include "Molecule.h"

// Canonical species string helpers
// Currently used for public DB (Ensembl) paths; extend as needed
inline const char* speciesToString(Species species)
{
    switch (species)
    {
    case Species::GENERIC: return "generic";
    case Species::C_ELEGANS: return "caenorhabditis_elegans";
    default:
        assert(false && "Unsupported species in speciesToString");
        return "unknown";
    }
}

inline std::wstring speciesToWString(Species species)
{
    const char* s = speciesToString(species);
    return std::wstring(s, s + std::char_traits<char>::length(s));
}


