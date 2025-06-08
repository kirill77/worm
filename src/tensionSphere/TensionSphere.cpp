#include "pch.h"
#include "TensionSphere.h"
#include <cmath>
#include <algorithm>
#include <cassert>

constexpr double PI = 3.14159265358979323846;

// TensionSphere implementation
TensionSphere::TensionSphere(uint32_t subdivisionLevel)
{
    // Create the base mesh with icosahedron and subdivisions
    m_pMesh = std::make_shared<ConnectedMesh>(1, subdivisionLevel);
}

void TensionSphere::makeTimeStep(double fDtSec)
{
    fDtSec = fDtSec;
}