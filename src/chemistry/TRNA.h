#pragma once

#include "StringDict.h"
#include <string>
#include <array>
#include <vector>

// Utility class for tRNA ID operations
// This class contains only static methods for working with tRNA molecule IDs
class TRNA
{
public:
    // Check if an ID represents a charged tRNA
    static bool isChargedTRNA(StringDict::ID id);
    
    // Convert an uncharged tRNA ID to its charged variant
    static StringDict::ID getChargedVariant(StringDict::ID unchargedID);
    
    // Get the anticodon sequence for a tRNA based on its StringDict ID
    static std::string getAnticodon(StringDict::ID tRNAId);
    
    // Get charged tRNA IDs that have a specific anticodon (for codon matching)
    static std::vector<StringDict::ID> getChargedTRNAsWithAnticodon(const std::string& anticodon);
    
    // Helper to convert codon to anticodon
    static std::string codonToAnticodon(const std::string& codon);
    
    // Get array of all uncharged tRNA IDs
    static const std::array<StringDict::ID, 25>& getUnchargedTRNAIds();
    
    // Test all TRNA functionality (called during initialization)
    static void runTests();

private:
    // Private constructor - this is a utility class with only static methods
    TRNA() = delete;
};
