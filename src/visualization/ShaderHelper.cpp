#include "pch.h"
#include "ShaderHelper.h"
#include <fstream>
#include <vector>
#include <d3dcompiler.h>

#pragma comment(lib, "d3dcompiler.lib")

ShaderHelper& ShaderHelper::getInstance()
{
    static ShaderHelper instance;
    return instance;
}

Microsoft::WRL::ComPtr<ID3DBlob> ShaderHelper::loadShader(
    const std::wstring& filePath,
    const std::string& entryPoint,
    const std::string& target,
    UINT compileFlags)
{
    // Create a unique key for the shader
    std::wstring shaderKey = filePath + L":" + std::wstring(entryPoint.begin(), entryPoint.end()) + L":" + std::wstring(target.begin(), target.end());
    
    // Check if already loaded
    auto it = m_shaderCache.find(shaderKey);
    if (it != m_shaderCache.end())
    {
        return it->second;
    }
    
    // Compile shader
    Microsoft::WRL::ComPtr<ID3DBlob> shaderBlob;
    Microsoft::WRL::ComPtr<ID3DBlob> errorBlob;
    
    HRESULT hr = D3DCompileFromFile(
        filePath.c_str(),
        nullptr,
        D3D_COMPILE_STANDARD_FILE_INCLUDE,
        entryPoint.c_str(),
        target.c_str(),
        compileFlags,
        0,
        &shaderBlob,
        &errorBlob);
    
    if (FAILED(hr))
    {
        if (errorBlob)
        {
            OutputDebugStringA(static_cast<const char*>(errorBlob->GetBufferPointer()));
        }
        ThrowIfFailed(hr);
    }
    
    // Cache the shader
    m_shaderCache[shaderKey] = shaderBlob;
    
    return shaderBlob;
}

Microsoft::WRL::ComPtr<ID3DBlob> ShaderHelper::loadCompiledShader(const std::wstring& filePath)
{
    // Check if already loaded
    auto it = m_shaderCache.find(filePath);
    if (it != m_shaderCache.end())
    {
        return it->second;
    }
    
    // Load the shader
    Microsoft::WRL::ComPtr<ID3DBlob> shaderBlob;
    ThrowIfFailed(D3DReadFileToBlob(filePath.c_str(), &shaderBlob));
    
    // Cache the shader
    m_shaderCache[filePath] = shaderBlob;
    
    return shaderBlob;
}

void ShaderHelper::clearCache()
{
    m_shaderCache.clear();
} 