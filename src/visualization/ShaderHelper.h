#pragma once

#include "pch.h"
#include "DirectXHelpers.h"
#include <string>
#include <unordered_map>

class ShaderHelper
{
public:
    // Singleton instance
    static ShaderHelper& getInstance();

    // Load a shader from file or cache
    Microsoft::WRL::ComPtr<ID3DBlob> loadShader(
        const std::wstring& filePath,
        const std::string& entryPoint,
        const std::string& target,
        UINT compileFlags = 0);
    
    // Load precompiled shader from file
    Microsoft::WRL::ComPtr<ID3DBlob> loadCompiledShader(const std::wstring& filePath);
    
    // Clear the shader cache
    void clearCache();

private:
    // Private constructor for singleton
    ShaderHelper() {}
    
    // Delete copy and move constructors and assign operators
    ShaderHelper(const ShaderHelper&) = delete;
    ShaderHelper& operator=(const ShaderHelper&) = delete;
    ShaderHelper(ShaderHelper&&) = delete;
    ShaderHelper& operator=(ShaderHelper&&) = delete;

    // Shader cache
    std::unordered_map<std::wstring, Microsoft::WRL::ComPtr<ID3DBlob>> m_shaderCache;
}; 