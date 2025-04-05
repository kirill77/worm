#include "pch.h"
#include "GPUStats.h"
#include "DirectXHelpers.h"
#include "GPUQueue.h"
#include <sstream>
#include <iomanip>

GPUStats::GPUStats(Microsoft::WRL::ComPtr<ID3D12Device> device, std::shared_ptr<GPUQueue> gpuQueue)
    : m_device(device)
    , m_gpuQueue(gpuQueue)
    , m_isCollecting(false)
{
    // Get timestamp frequency
    m_timestampFrequency = 0;
    m_gpuQueue->getQueue()->GetTimestampFrequency(&m_timestampFrequency);
    
    initializeQueries();
}

GPUStats::~GPUStats()
{
    // Ensure we're not collecting stats when destroyed
    if (m_isCollecting)
    {
        end();
    }
}

void GPUStats::initializeQueries()
{
    // Create query heap for pipeline statistics
    D3D12_QUERY_HEAP_DESC pipelineStatsQueryHeapDesc = {};
    pipelineStatsQueryHeapDesc.Type = D3D12_QUERY_HEAP_TYPE_PIPELINE_STATISTICS;
    pipelineStatsQueryHeapDesc.Count = 2; // Begin and end queries
    pipelineStatsQueryHeapDesc.NodeMask = 0;
    
    ThrowIfFailed(m_device->CreateQueryHeap(&pipelineStatsQueryHeapDesc, IID_PPV_ARGS(&m_queryHeap)));
    
    // Create query buffer
    D3D12_RESOURCE_DESC queryBufferDesc = {};
    queryBufferDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    queryBufferDesc.Width = sizeof(D3D12_QUERY_DATA_PIPELINE_STATISTICS) * 2; // Space for begin and end queries
    queryBufferDesc.Height = 1;
    queryBufferDesc.DepthOrArraySize = 1;
    queryBufferDesc.MipLevels = 1;
    queryBufferDesc.Format = DXGI_FORMAT_UNKNOWN;
    queryBufferDesc.SampleDesc.Count = 1;
    queryBufferDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
    queryBufferDesc.Flags = D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS;
    
    CD3DX12_HEAP_PROPERTIES defaultHeapProperties(D3D12_HEAP_TYPE_DEFAULT);
    ThrowIfFailed(m_device->CreateCommittedResource(
        &defaultHeapProperties,
        D3D12_HEAP_FLAG_NONE,
        &queryBufferDesc,
        D3D12_RESOURCE_STATE_COPY_DEST,
        nullptr,
        IID_PPV_ARGS(&m_queryBuffer)));
    
    // Create readback buffer
    CD3DX12_HEAP_PROPERTIES readbackHeapProperties(D3D12_HEAP_TYPE_READBACK);
    ThrowIfFailed(m_device->CreateCommittedResource(
        &readbackHeapProperties,
        D3D12_HEAP_FLAG_NONE,
        &queryBufferDesc,
        D3D12_RESOURCE_STATE_COPY_DEST,
        nullptr,
        IID_PPV_ARGS(&m_queryReadbackBuffer)));
    
    // Initialize query indices
    m_queryIndices.pipelineStatsBegin = 0;
    m_queryIndices.pipelineStatsEnd = 1;
    m_queryIndices.timestampBegin = 0;
    m_queryIndices.timestampEnd = 1;
}

void GPUStats::begin()
{
    if (m_isCollecting)
    {
        return;
    }
    
    // Get command list from queue
    auto commandList = m_gpuQueue->beginRecording();
    
    // Begin pipeline statistics query
    commandList->BeginQuery(m_queryHeap.Get(), D3D12_QUERY_TYPE_PIPELINE_STATISTICS, m_queryIndices.pipelineStatsBegin);
    
    // Insert timestamp
    commandList->EndQuery(m_queryHeap.Get(), D3D12_QUERY_TYPE_TIMESTAMP, m_queryIndices.timestampBegin);
    
    // Close and execute command list
    commandList->Close();
    m_gpuQueue->execute(commandList);
    
    m_isCollecting = true;
}

void GPUStats::end()
{
    if (!m_isCollecting)
    {
        return;
    }
    
    // Get command list from queue
    auto commandList = m_gpuQueue->beginRecording();
    
    // End pipeline statistics query
    commandList->EndQuery(m_queryHeap.Get(), D3D12_QUERY_TYPE_PIPELINE_STATISTICS, m_queryIndices.pipelineStatsEnd);
    
    // Insert timestamp
    commandList->EndQuery(m_queryHeap.Get(), D3D12_QUERY_TYPE_TIMESTAMP, m_queryIndices.timestampEnd);
    
    // Resolve query data
    commandList->ResolveQueryData(
        m_queryHeap.Get(),
        D3D12_QUERY_TYPE_PIPELINE_STATISTICS,
        m_queryIndices.pipelineStatsBegin,
        2, // Number of queries to resolve
        m_queryBuffer.Get(),
        0);
    
    // Create resource barrier
    auto barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        m_queryBuffer.Get(),
        D3D12_RESOURCE_STATE_COPY_DEST,
        D3D12_RESOURCE_STATE_COPY_SOURCE);
    
    // Copy query data to readback buffer
    commandList->ResourceBarrier(1, &barrier);
    
    commandList->CopyResource(m_queryReadbackBuffer.Get(), m_queryBuffer.Get());
    
    // Close and execute command list
    commandList->Close();
    m_gpuQueue->execute(commandList);
    
    // Read back query data
    readQueryData();
    
    m_isCollecting = false;
}

void GPUStats::readQueryData()
{
    // Map the readback buffer
    D3D12_QUERY_DATA_PIPELINE_STATISTICS* pData = nullptr;
    CD3DX12_RANGE readRange(0, sizeof(D3D12_QUERY_DATA_PIPELINE_STATISTICS));
    ThrowIfFailed(m_queryReadbackBuffer->Map(0, &readRange, reinterpret_cast<void**>(&pData)));
    
    // Copy the data
    m_pipelineStats = *pData;
    
    // Unmap the buffer
    m_queryReadbackBuffer->Unmap(0, nullptr);
    
    // Read timestamp data
    uint64_t* pTimestampData = nullptr;
    ThrowIfFailed(m_queryReadbackBuffer->Map(0, &readRange, reinterpret_cast<void**>(&pTimestampData)));
    
    m_timestampBegin = pTimestampData[0];
    m_timestampEnd = pTimestampData[1];
    
    m_queryReadbackBuffer->Unmap(0, nullptr);
}

std::string GPUStats::getStats() const
{
    std::stringstream ss;
    
    // Calculate time in milliseconds
    double timeMs = 0.0;
    if (m_timestampFrequency > 0)
    {
        timeMs = static_cast<double>(m_timestampEnd - m_timestampBegin) / static_cast<double>(m_timestampFrequency) * 1000.0;
    }
    
    ss << std::fixed << std::setprecision(2);
    ss << "GPU Statistics:\n";
    ss << "Time: " << timeMs << " ms\n";
    ss << "Vertex Shader Invocations: " << m_pipelineStats.IAVertices << "\n";
    ss << "Input Assembler Primitives: " << m_pipelineStats.IAPrimitives << "\n";
    ss << "Vertex Shader Invocations: " << m_pipelineStats.VSInvocations << "\n";
    ss << "Geometry Shader Invocations: " << m_pipelineStats.GSInvocations << "\n";
    ss << "Geometry Shader Primitives: " << m_pipelineStats.GSPrimitives << "\n";
    ss << "Clipping Stage Primitives: " << m_pipelineStats.CInvocations << "\n";
    ss << "Clipping Stage Primitives: " << m_pipelineStats.CPrimitives << "\n";
    ss << "Pixel Shader Invocations: " << m_pipelineStats.PSInvocations << "\n";
    ss << "Hull Shader Invocations: " << m_pipelineStats.HSInvocations << "\n";
    ss << "Domain Shader Invocations: " << m_pipelineStats.DSInvocations << "\n";
    ss << "Compute Shader Invocations: " << m_pipelineStats.CSInvocations << "\n";
    
    return ss.str();
} 