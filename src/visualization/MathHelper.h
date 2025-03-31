#pragma once

// Custom min/max functions to avoid conflicts with Windows.h macros
// These will be used in our visualization code instead of relying on std::min/max

template <typename T>
inline T mathMin(T a, T b) {
    return (a < b) ? a : b;
}

template <typename T>
inline T mathMax(T a, T b) {
    return (a > b) ? a : b;
} 