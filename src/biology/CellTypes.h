#pragma once

enum class CellType
{
    Zygote,         // First cell (formerly P0)
    AB,             // Anterior blastomere
    Germline1,      // First germline cell (formerly P1)
    EMS,            // Endoderm/Mesoderm precursor
    Germline2,      // Second germline cell (formerly P2)
    Germline3,      // Third germline cell (formerly P3)
    PrimordialGerm, // Primordial germ cell (formerly P4)
    Somatic         // Other somatic cells
}; 