// pch.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing features.
// However, files listed here are ALL re-compiled if any one of them is updated between builds.
// Do not add files here that you will be updating frequently as this negates the performance advantage.

#ifndef PCH_H
#define PCH_H

// Define feature level for required D3D12 features
// D3D12_SDK_VERSION is defined in d3d12.h - don't redefine it here
#define D3D12_FEATURE_LEVEL D3D12_FEATURE_LEVEL_12_0

// Prevent Windows.h min/max macro conflicts
#define NOMINMAX

// add headers that you want to pre-compile here
#include "framework.h"

// Windows and DirectX headers
#include <windows.h>
#include <windowsx.h>
#include <wrl.h>
#include <d3d12.h>
#include <dxgi1_6.h>
#include <d3dcompiler.h>
#include <DirectXMath.h>
#include <DirectXColors.h>
#include "CD3DX12.h"
#include "MathHelper.h"

// STL headers
#include <memory>
#include <vector>
#include <string>
#include <array>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <exception>

// C++20 headers - keep these together so they can be conditionally included
#if _MSVC_LANG >= 202002L  // C++20 or later
#include <numbers>
#include <concepts>
#endif

// Link necessary libraries
#pragma comment(lib, "d3d12.lib")
#pragma comment(lib, "dxgi.lib")
#pragma comment(lib, "d3dcompiler.lib")

#endif //PCH_H
