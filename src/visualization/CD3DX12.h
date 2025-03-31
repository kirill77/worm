//*********************************************************
//
// Copyright (c) Microsoft. All rights reserved.
// This code is licensed under the MIT License (MIT).
// THIS CODE IS PROVIDED *AS IS* WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING ANY
// IMPLIED WARRANTIES OF FITNESS FOR A PARTICULAR
// PURPOSE, MERCHANTABILITY, OR NON-INFRINGEMENT.
//
//*********************************************************

// This file is a simplified version of the CD3DX12 helpers from the DirectX samples

#pragma once

#include "pch.h"

// Fallback to Windows SDK definitions if D3D12 helper library is not included
// This is a simplified version for our needs - not a full replacement of d3dx12.h

//================================================================================================
// D3DX12 Resource Barrier Helper
//================================================================================================

struct CD3DX12_RESOURCE_BARRIER : public D3D12_RESOURCE_BARRIER
{
    CD3DX12_RESOURCE_BARRIER() = default;
    
    static inline CD3DX12_RESOURCE_BARRIER Transition(
        _In_ ID3D12Resource* pResource,
        D3D12_RESOURCE_STATES stateBefore,
        D3D12_RESOURCE_STATES stateAfter,
        UINT subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES)
    {
        CD3DX12_RESOURCE_BARRIER result;
        ZeroMemory(&result, sizeof(result));
        D3D12_RESOURCE_BARRIER &barrier = result;
        barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
        barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
        barrier.Transition.pResource = pResource;
        barrier.Transition.StateBefore = stateBefore;
        barrier.Transition.StateAfter = stateAfter;
        barrier.Transition.Subresource = subresource;
        return result;
    }
};

//================================================================================================
// Descriptor Range Helper
//================================================================================================

struct CD3DX12_DESCRIPTOR_RANGE : public D3D12_DESCRIPTOR_RANGE
{
    CD3DX12_DESCRIPTOR_RANGE() = default;
    
    inline CD3DX12_DESCRIPTOR_RANGE(
        D3D12_DESCRIPTOR_RANGE_TYPE rangeType,
        UINT numDescriptors,
        UINT baseShaderRegister,
        UINT registerSpace = 0,
        UINT offsetInDescriptorsFromTableStart = D3D12_DESCRIPTOR_RANGE_OFFSET_APPEND)
    {
        Init(rangeType, numDescriptors, baseShaderRegister, registerSpace, offsetInDescriptorsFromTableStart);
    }
    
    inline void Init(
        D3D12_DESCRIPTOR_RANGE_TYPE rangeType,
        UINT numDescriptors,
        UINT baseShaderRegister,
        UINT registerSpace = 0,
        UINT offsetInDescriptorsFromTableStart = D3D12_DESCRIPTOR_RANGE_OFFSET_APPEND)
    {
        RangeType = rangeType;
        NumDescriptors = numDescriptors;
        BaseShaderRegister = baseShaderRegister;
        RegisterSpace = registerSpace;
        OffsetInDescriptorsFromTableStart = offsetInDescriptorsFromTableStart;
    }
};

//================================================================================================
// Root Parameter Helper
//================================================================================================

struct CD3DX12_ROOT_PARAMETER : public D3D12_ROOT_PARAMETER
{
    CD3DX12_ROOT_PARAMETER() = default;
    
    inline void InitAsDescriptorTable(
        UINT numDescriptorRanges,
        const D3D12_DESCRIPTOR_RANGE* pDescriptorRanges,
        D3D12_SHADER_VISIBILITY visibility = D3D12_SHADER_VISIBILITY_ALL)
    {
        ParameterType = D3D12_ROOT_PARAMETER_TYPE_DESCRIPTOR_TABLE;
        DescriptorTable.NumDescriptorRanges = numDescriptorRanges;
        DescriptorTable.pDescriptorRanges = pDescriptorRanges;
        ShaderVisibility = visibility;
    }
    
    inline void InitAsConstants(
        UINT num32BitValues,
        UINT shaderRegister,
        UINT registerSpace = 0,
        D3D12_SHADER_VISIBILITY visibility = D3D12_SHADER_VISIBILITY_ALL)
    {
        ParameterType = D3D12_ROOT_PARAMETER_TYPE_32BIT_CONSTANTS;
        Constants.Num32BitValues = num32BitValues;
        Constants.ShaderRegister = shaderRegister;
        Constants.RegisterSpace = registerSpace;
        ShaderVisibility = visibility;
    }
    
    inline void InitAsConstantBufferView(
        UINT shaderRegister,
        UINT registerSpace = 0,
        D3D12_SHADER_VISIBILITY visibility = D3D12_SHADER_VISIBILITY_ALL)
    {
        ParameterType = D3D12_ROOT_PARAMETER_TYPE_CBV;
        Descriptor.ShaderRegister = shaderRegister;
        Descriptor.RegisterSpace = registerSpace;
        ShaderVisibility = visibility;
    }
    
    inline void InitAsShaderResourceView(
        UINT shaderRegister,
        UINT registerSpace = 0,
        D3D12_SHADER_VISIBILITY visibility = D3D12_SHADER_VISIBILITY_ALL)
    {
        ParameterType = D3D12_ROOT_PARAMETER_TYPE_SRV;
        Descriptor.ShaderRegister = shaderRegister;
        Descriptor.RegisterSpace = registerSpace;
        ShaderVisibility = visibility;
    }
    
    inline void InitAsUnorderedAccessView(
        UINT shaderRegister,
        UINT registerSpace = 0,
        D3D12_SHADER_VISIBILITY visibility = D3D12_SHADER_VISIBILITY_ALL)
    {
        ParameterType = D3D12_ROOT_PARAMETER_TYPE_UAV;
        Descriptor.ShaderRegister = shaderRegister;
        Descriptor.RegisterSpace = registerSpace;
        ShaderVisibility = visibility;
    }
};

//================================================================================================
// Root Signature Helper
//================================================================================================

struct CD3DX12_ROOT_SIGNATURE_DESC : public D3D12_ROOT_SIGNATURE_DESC
{
    CD3DX12_ROOT_SIGNATURE_DESC() = default;
    
    inline CD3DX12_ROOT_SIGNATURE_DESC(
        UINT numParameters,
        const D3D12_ROOT_PARAMETER* _pParameters,
        UINT numStaticSamplers = 0,
        const D3D12_STATIC_SAMPLER_DESC* _pStaticSamplers = nullptr,
        D3D12_ROOT_SIGNATURE_FLAGS flags = D3D12_ROOT_SIGNATURE_FLAG_NONE)
    {
        Init(numParameters, _pParameters, numStaticSamplers, _pStaticSamplers, flags);
    }
    
    inline void Init(
        UINT numParameters,
        const D3D12_ROOT_PARAMETER* _pParameters,
        UINT numStaticSamplers = 0,
        const D3D12_STATIC_SAMPLER_DESC* _pStaticSamplers = nullptr,
        D3D12_ROOT_SIGNATURE_FLAGS flags = D3D12_ROOT_SIGNATURE_FLAG_NONE)
    {
        NumParameters = numParameters;
        pParameters = _pParameters;
        NumStaticSamplers = numStaticSamplers;
        pStaticSamplers = _pStaticSamplers;
        Flags = flags;
    }
};

//================================================================================================
// Pipeline State Stream Helper
//================================================================================================

struct CD3DX12_PIPELINE_STATE_STREAM_ROOT_SIGNATURE
{
    CD3DX12_PIPELINE_STATE_STREAM_ROOT_SIGNATURE() = default;
    explicit CD3DX12_PIPELINE_STATE_STREAM_ROOT_SIGNATURE(ID3D12RootSignature* pRootSignature) : pRootSignature(pRootSignature) { }
    ID3D12RootSignature* pRootSignature;
};

struct CD3DX12_PIPELINE_STATE_STREAM_INPUT_LAYOUT
{
    CD3DX12_PIPELINE_STATE_STREAM_INPUT_LAYOUT() = default;
    CD3DX12_PIPELINE_STATE_STREAM_INPUT_LAYOUT(const D3D12_INPUT_LAYOUT_DESC& inputLayoutDesc) : InputLayout(inputLayoutDesc) { }
    D3D12_INPUT_LAYOUT_DESC InputLayout;
};

struct CD3DX12_PIPELINE_STATE_STREAM_PRIMITIVE_TOPOLOGY
{
    CD3DX12_PIPELINE_STATE_STREAM_PRIMITIVE_TOPOLOGY() = default;
    explicit CD3DX12_PIPELINE_STATE_STREAM_PRIMITIVE_TOPOLOGY(D3D12_PRIMITIVE_TOPOLOGY_TYPE primitiveTopologyType) : PrimitiveTopologyType(primitiveTopologyType) { }
    D3D12_PRIMITIVE_TOPOLOGY_TYPE PrimitiveTopologyType;
};

struct CD3DX12_PIPELINE_STATE_STREAM_VS
{
    CD3DX12_PIPELINE_STATE_STREAM_VS() = default;
    explicit CD3DX12_PIPELINE_STATE_STREAM_VS(const D3D12_SHADER_BYTECODE& vs) : VS(vs) { }
    D3D12_SHADER_BYTECODE VS;
};

struct CD3DX12_PIPELINE_STATE_STREAM_PS
{
    CD3DX12_PIPELINE_STATE_STREAM_PS() = default;
    explicit CD3DX12_PIPELINE_STATE_STREAM_PS(const D3D12_SHADER_BYTECODE& ps) : PS(ps) { }
    D3D12_SHADER_BYTECODE PS;
};

struct CD3DX12_PIPELINE_STATE_STREAM_DEPTH_STENCIL_FORMAT
{
    CD3DX12_PIPELINE_STATE_STREAM_DEPTH_STENCIL_FORMAT() = default;
    explicit CD3DX12_PIPELINE_STATE_STREAM_DEPTH_STENCIL_FORMAT(DXGI_FORMAT dsvFormat) : DSVFormat(dsvFormat) { }
    DXGI_FORMAT DSVFormat;
};

struct CD3DX12_PIPELINE_STATE_STREAM_RENDER_TARGET_FORMATS
{
    CD3DX12_PIPELINE_STATE_STREAM_RENDER_TARGET_FORMATS() = default;
    explicit CD3DX12_PIPELINE_STATE_STREAM_RENDER_TARGET_FORMATS(const D3D12_RT_FORMAT_ARRAY& formats) : RTFormats(formats) { }
    D3D12_RT_FORMAT_ARRAY RTFormats;
};

//================================================================================================
// Shader Bytecode Helper
//================================================================================================

struct CD3DX12_SHADER_BYTECODE : public D3D12_SHADER_BYTECODE
{
    CD3DX12_SHADER_BYTECODE() = default;
    
    explicit CD3DX12_SHADER_BYTECODE(const D3D12_SHADER_BYTECODE& o) :
        D3D12_SHADER_BYTECODE(o)
    {}
    
    CD3DX12_SHADER_BYTECODE(
        _In_reads_bytes_(bytecodeLength) const void* _pShaderBytecode,
        SIZE_T bytecodeLength)
    {
        pShaderBytecode = _pShaderBytecode;
        BytecodeLength = bytecodeLength;
    }
};

//================================================================================================
// Heap Properties Helper
//================================================================================================

struct CD3DX12_HEAP_PROPERTIES : public D3D12_HEAP_PROPERTIES
{
    CD3DX12_HEAP_PROPERTIES() = default;
    
    explicit CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE type)
    {
        Type = type;
        CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
        MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
        CreationNodeMask = 1;
        VisibleNodeMask = 1;
    }
};

//================================================================================================
// Resource Desc Helper
//================================================================================================

struct CD3DX12_RESOURCE_DESC : public D3D12_RESOURCE_DESC
{
    CD3DX12_RESOURCE_DESC()
    {
        Dimension = D3D12_RESOURCE_DIMENSION_UNKNOWN;
        Alignment = 0;
        Width = 0;
        Height = 1;
        DepthOrArraySize = 1;
        MipLevels = 1;
        Format = DXGI_FORMAT_UNKNOWN;
        SampleDesc.Count = 1;
        SampleDesc.Quality = 0;
        Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN;
        Flags = D3D12_RESOURCE_FLAG_NONE;
    }
    
    static inline CD3DX12_RESOURCE_DESC Buffer(
        UINT64 width,
        D3D12_RESOURCE_FLAGS flags = D3D12_RESOURCE_FLAG_NONE,
        UINT64 alignment = 0)
    {
        CD3DX12_RESOURCE_DESC result;
        result.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
        result.Alignment = alignment;
        result.Width = width;
        result.Height = 1;
        result.DepthOrArraySize = 1;
        result.MipLevels = 1;
        result.Format = DXGI_FORMAT_UNKNOWN;
        result.SampleDesc.Count = 1;
        result.SampleDesc.Quality = 0;
        result.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
        result.Flags = flags;
        return result;
    }
};

//================================================================================================
// CPU and GPU Descriptor Handle Helpers
//================================================================================================

struct CD3DX12_CPU_DESCRIPTOR_HANDLE : public D3D12_CPU_DESCRIPTOR_HANDLE
{
    CD3DX12_CPU_DESCRIPTOR_HANDLE() = default;
    explicit CD3DX12_CPU_DESCRIPTOR_HANDLE(D3D12_CPU_DESCRIPTOR_HANDLE handle) : D3D12_CPU_DESCRIPTOR_HANDLE(handle) {}
    
    CD3DX12_CPU_DESCRIPTOR_HANDLE(
        D3D12_CPU_DESCRIPTOR_HANDLE handle, 
        INT offsetScaledByIncrementSize)
    {
        InitOffsetted(handle, offsetScaledByIncrementSize);
    }
    
    CD3DX12_CPU_DESCRIPTOR_HANDLE(
        D3D12_CPU_DESCRIPTOR_HANDLE handle,
        INT offsetInDescriptors,
        UINT descriptorIncrementSize)
    {
        InitOffsetted(handle, offsetInDescriptors, descriptorIncrementSize);
    }
    
    void InitOffsetted(
        D3D12_CPU_DESCRIPTOR_HANDLE handle, 
        INT offsetScaledByIncrementSize)
    {
        ptr = handle.ptr + offsetScaledByIncrementSize;
    }
    
    void InitOffsetted(
        D3D12_CPU_DESCRIPTOR_HANDLE handle,
        INT offsetInDescriptors,
        UINT descriptorIncrementSize)
    {
        ptr = handle.ptr + offsetInDescriptors * descriptorIncrementSize;
    }
};

struct CD3DX12_GPU_DESCRIPTOR_HANDLE : public D3D12_GPU_DESCRIPTOR_HANDLE
{
    CD3DX12_GPU_DESCRIPTOR_HANDLE() = default;
    explicit CD3DX12_GPU_DESCRIPTOR_HANDLE(D3D12_GPU_DESCRIPTOR_HANDLE handle) : D3D12_GPU_DESCRIPTOR_HANDLE(handle) {}
    
    CD3DX12_GPU_DESCRIPTOR_HANDLE(
        D3D12_GPU_DESCRIPTOR_HANDLE handle,
        INT offsetInDescriptors,
        UINT descriptorIncrementSize)
    {
        InitOffsetted(handle, offsetInDescriptors, descriptorIncrementSize);
    }
    
    void InitOffsetted(
        D3D12_GPU_DESCRIPTOR_HANDLE handle, 
        INT offsetInDescriptors,
        UINT descriptorIncrementSize)
    {
        ptr = handle.ptr + offsetInDescriptors * descriptorIncrementSize;
    }
};

//================================================================================================
// Range Helper
//================================================================================================

struct CD3DX12_RANGE : public D3D12_RANGE
{
    CD3DX12_RANGE() = default;
    CD3DX12_RANGE(SIZE_T begin, SIZE_T end) 
    { 
        Begin = begin; 
        End = end; 
    }
};

//================================================================================================
// UpdateSubresources Helper
//================================================================================================

inline UINT64 UpdateSubresources(
    _In_ ID3D12GraphicsCommandList* pCmdList,
    _In_ ID3D12Resource* pDestinationResource,
    _In_ ID3D12Resource* pIntermediate,
    _In_range_(0,D3D12_REQ_SUBRESOURCES) UINT64 FirstSubresource,
    _In_range_(0,D3D12_REQ_SUBRESOURCES-FirstSubresource) UINT64 NumSubresources,
    UINT64 RequiredSize,
    _In_reads_(NumSubresources) const D3D12_PLACED_SUBRESOURCE_FOOTPRINT* pLayouts,
    _In_reads_(NumSubresources) const UINT* pNumRows,
    _In_reads_(NumSubresources) const UINT64* pRowSizesInBytes,
    _In_reads_(NumSubresources) const D3D12_SUBRESOURCE_DATA* pSrcData)
{
    // Minor validation
    if (pCmdList == nullptr || pDestinationResource == nullptr || pIntermediate == nullptr || pLayouts == nullptr || pNumRows == nullptr || pRowSizesInBytes == nullptr || pSrcData == nullptr)
        return 0;
        
    UINT64 MemToAlloc = static_cast<UINT64>(sizeof(D3D12_PLACED_SUBRESOURCE_FOOTPRINT) + sizeof(UINT) + sizeof(UINT64)) * NumSubresources;
    if (MemToAlloc > SIZE_MAX)
        return 0;
        
    // Update subresources
    UINT64 IntermediateOffset = 0;
    UINT64 IntermediateSize = 0;
    UINT64 TotalSize = 0;
    
    for (UINT i = 0; i < FirstSubresource; ++i)
    {
        TotalSize += pRowSizesInBytes[i] * pNumRows[i];
    }
    
    // Copy data to intermediate resource
    BYTE* pData;
    pIntermediate->Map(0, nullptr, reinterpret_cast<void**>(&pData));
    
    for (UINT i = 0; i < NumSubresources; ++i)
    {
        UINT64 RowSizeInBytes = pRowSizesInBytes[i];
        UINT RowCount = pNumRows[i];
        
        if (RowSizeInBytes > SIZE_MAX || RowCount > SIZE_MAX)
            return 0;
            
        const D3D12_SUBRESOURCE_DATA& srcData = pSrcData[i];
        const BYTE* pSrcRow = static_cast<const BYTE*>(srcData.pData);
        BYTE* pDestRow = pData + pLayouts[i].Offset;
        
        for (UINT j = 0; j < RowCount; ++j)
        {
            memcpy(pDestRow, pSrcRow, static_cast<SIZE_T>(RowSizeInBytes));
            pSrcRow += srcData.RowPitch;
            pDestRow += pLayouts[i].Footprint.RowPitch;
        }
        
        TotalSize += RowSizeInBytes * RowCount;
    }
    
    pIntermediate->Unmap(0, nullptr);
    
    // Schedule copy from intermediate to destination
    pCmdList->CopyResource(pDestinationResource, pIntermediate);
    
    return RequiredSize;
}

inline UINT64 UpdateSubresources(
    _In_ ID3D12GraphicsCommandList* pCmdList,
    _In_ ID3D12Resource* pDestinationResource,
    _In_ ID3D12Resource* pIntermediate,
    UINT64 IntermediateOffset,
    _In_range_(0,D3D12_REQ_SUBRESOURCES) UINT FirstSubresource,
    _In_range_(0,D3D12_REQ_SUBRESOURCES-FirstSubresource) UINT NumSubresources,
    _In_reads_(NumSubresources) const D3D12_SUBRESOURCE_DATA* pSrcData)
{
    UINT64 RequiredSize = 0;
    UINT64 MemToAlloc = static_cast<UINT64>(sizeof(D3D12_PLACED_SUBRESOURCE_FOOTPRINT) + sizeof(UINT) + sizeof(UINT64)) * NumSubresources;
    if (MemToAlloc > SIZE_MAX)
        return 0;
    
    void* pMem = HeapAlloc(GetProcessHeap(), 0, static_cast<SIZE_T>(MemToAlloc));
    if (pMem == nullptr)
        return 0;
    
    D3D12_PLACED_SUBRESOURCE_FOOTPRINT* pLayouts = static_cast<D3D12_PLACED_SUBRESOURCE_FOOTPRINT*>(pMem);
    UINT* pNumRows = reinterpret_cast<UINT*>(pLayouts + NumSubresources);
    UINT64* pRowSizesInBytes = reinterpret_cast<UINT64*>(pNumRows + NumSubresources);
    
    D3D12_RESOURCE_DESC Desc = pDestinationResource->GetDesc();
    ID3D12Device* pDevice;
    pDestinationResource->GetDevice(IID_PPV_ARGS(&pDevice));
    pDevice->GetCopyableFootprints(&Desc, FirstSubresource, NumSubresources, IntermediateOffset, pLayouts, pNumRows, pRowSizesInBytes, &RequiredSize);
    pDevice->Release();
    
    UINT64 Result = UpdateSubresources(pCmdList, pDestinationResource, pIntermediate, FirstSubresource, NumSubresources, RequiredSize, pLayouts, pNumRows, pRowSizesInBytes, pSrcData);
    HeapFree(GetProcessHeap(), 0, pMem);
    return Result;
} 