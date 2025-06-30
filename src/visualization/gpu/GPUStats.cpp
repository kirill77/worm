#include "pch.h"
#include "GPUStats.h"
#include "DirectXHelpers.h"
#include "GPUQueue.h"
#include <sstream>
#include <iomanip>

GPUStats::GPUStats(Microsoft::WRL::ComPtr<ID3D12Device> device)
    : m_device(device)
    , m_isCollecting(false)
{
    // Get timestamp frequency from the device
    m_timestampFrequency = 0;
    
    // Create a temporary command queue to get the timestamp frequency
    D3D12_COMMAND_QUEUE_DESC queueDesc = {};
    queueDesc.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;
    queueDesc.Flags = D3D12_COMMAND_QUEUE_FLAG_NONE;
    
    Microsoft::WRL::ComPtr<ID3D12CommandQueue> tempQueue;
    ThrowIfFailed(m_device->CreateCommandQueue(&queueDesc, IID_PPV_ARGS(&tempQueue)));
    tempQueue->GetTimestampFrequency(&m_timestampFrequency);
    
    initializeQueries();
}

GPUStats::~GPUStats()
{
    // Ensure we're not collecting stats when destroyed
    if (m_isCollecting)
    {
        // We can't call end() here because we don't have a command list
        m_isCollecting = false;
    }
}

void GPUStats::initializeQueries()
{
    // Create query heap for pipeline statistics
    D3D12_QUERY_HEAP_DESC pipelineStatsQueryHeapDesc = {};
    pipelineStatsQueryHeapDesc.Type = D3D12_QUERY_HEAP_TYPE_PIPELINE_STATISTICS;
    pipelineStatsQueryHeapDesc.Count = 2; // Begin and end queries
    pipelineStatsQueryHeapDesc.NodeMask = 0;
    
    ThrowIfFailed(m_device->CreateQueryHeap(&pipelineStatsQueryHeapDesc, IID_PPV_ARGS(&m_pipelineStatsQueryHeap)));
    
    // Create query heap for timestamps
    D3D12_QUERY_HEAP_DESC timestampQueryHeapDesc = {};
    timestampQueryHeapDesc.Type = D3D12_QUERY_HEAP_TYPE_TIMESTAMP;
    timestampQueryHeapDesc.Count = 2; // Begin and end queries
    timestampQueryHeapDesc.NodeMask = 0;
    
    ThrowIfFailed(m_device->CreateQueryHeap(&timestampQueryHeapDesc, IID_PPV_ARGS(&m_timestampQueryHeap)));
    
    // Create pipeline statistics query buffer
    D3D12_RESOURCE_DESC pipelineStatsBufferDesc = {};
    pipelineStatsBufferDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    pipelineStatsBufferDesc.Width = sizeof(D3D12_QUERY_DATA_PIPELINE_STATISTICS) * 2; // Space for begin and end queries
    pipelineStatsBufferDesc.Height = 1;
    pipelineStatsBufferDesc.DepthOrArraySize = 1;
    pipelineStatsBufferDesc.MipLevels = 1;
    pipelineStatsBufferDesc.Format = DXGI_FORMAT_UNKNOWN;
    pipelineStatsBufferDesc.SampleDesc.Count = 1;
    pipelineStatsBufferDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
    pipelineStatsBufferDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
    
    CD3DX12_HEAP_PROPERTIES defaultHeapProperties(D3D12_HEAP_TYPE_DEFAULT);
    ThrowIfFailed(m_device->CreateCommittedResource(
        &defaultHeapProperties,
        D3D12_HEAP_FLAG_NONE,
        &pipelineStatsBufferDesc,
        D3D12_RESOURCE_STATE_COPY_DEST,
        nullptr,
        IID_PPV_ARGS(&m_pipelineStatsBuffer)));
    
    // Create pipeline statistics readback buffer
    CD3DX12_HEAP_PROPERTIES readbackHeapProperties(D3D12_HEAP_TYPE_READBACK);
    ThrowIfFailed(m_device->CreateCommittedResource(
        &readbackHeapProperties,
        D3D12_HEAP_FLAG_NONE,
        &pipelineStatsBufferDesc,
        D3D12_RESOURCE_STATE_COPY_DEST,
        nullptr,
        IID_PPV_ARGS(&m_pipelineStatsReadbackBuffer)));
    
    // Create timestamp query buffer
    D3D12_RESOURCE_DESC timestampBufferDesc = {};
    timestampBufferDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    timestampBufferDesc.Width = sizeof(uint64_t) * 2; // Space for begin and end timestamps
    timestampBufferDesc.Height = 1;
    timestampBufferDesc.DepthOrArraySize = 1;
    timestampBufferDesc.MipLevels = 1;
    timestampBufferDesc.Format = DXGI_FORMAT_UNKNOWN;
    timestampBufferDesc.SampleDesc.Count = 1;
    timestampBufferDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
    timestampBufferDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
    
    ThrowIfFailed(m_device->CreateCommittedResource(
        &defaultHeapProperties,
        D3D12_HEAP_FLAG_NONE,
        &timestampBufferDesc,
        D3D12_RESOURCE_STATE_COPY_DEST,
        nullptr,
        IID_PPV_ARGS(&m_timestampBuffer)));
    
    // Create timestamp readback buffer
    ThrowIfFailed(m_device->CreateCommittedResource(
        &readbackHeapProperties,
        D3D12_HEAP_FLAG_NONE,
        &timestampBufferDesc,
        D3D12_RESOURCE_STATE_COPY_DEST,
        nullptr,
        IID_PPV_ARGS(&m_timestampReadbackBuffer)));
    
    // Initialize query indices
    m_queryIndices.pipelineStatsBegin = 0;
    m_queryIndices.pipelineStatsEnd = 0;
    m_queryIndices.timestampBegin = 0;
    m_queryIndices.timestampEnd = 1;
}

void GPUStats::begin(ID3D12GraphicsCommandList& cmdList)
{
    if (m_isCollecting)
    {
        return;
    }
    
    // Begin pipeline statistics query
    cmdList.BeginQuery(m_pipelineStatsQueryHeap.Get(), D3D12_QUERY_TYPE_PIPELINE_STATISTICS, m_queryIndices.pipelineStatsBegin);
    
    // Insert timestamp (no BeginQuery needed for timestamps)
    cmdList.EndQuery(m_timestampQueryHeap.Get(), D3D12_QUERY_TYPE_TIMESTAMP, m_queryIndices.timestampBegin);
    
    m_isCollecting = true;
}

void GPUStats::end(ID3D12GraphicsCommandList& cmdList)
{
    if (!m_isCollecting)
    {
        return;
    }
    
    // End pipeline statistics query
    cmdList.EndQuery(m_pipelineStatsQueryHeap.Get(), D3D12_QUERY_TYPE_PIPELINE_STATISTICS, m_queryIndices.pipelineStatsEnd);
    
    // Insert timestamp (no BeginQuery needed for timestamps)
    cmdList.EndQuery(m_timestampQueryHeap.Get(), D3D12_QUERY_TYPE_TIMESTAMP, m_queryIndices.timestampEnd);
    
    // Resolve pipeline statistics query data
    cmdList.ResolveQueryData(
        m_pipelineStatsQueryHeap.Get(),
        D3D12_QUERY_TYPE_PIPELINE_STATISTICS,
        m_queryIndices.pipelineStatsBegin,
        1, // Number of queries to resolve
        m_pipelineStatsBuffer.Get(),
        0);
    
    // Resolve timestamp query data
    cmdList.ResolveQueryData(
        m_timestampQueryHeap.Get(),
        D3D12_QUERY_TYPE_TIMESTAMP,
        m_queryIndices.timestampBegin,
        2, // Number of queries to resolve
        m_timestampBuffer.Get(),
        0);
    
    // Create resource barriers
    auto pipelineStatsBarrier = CD3DX12_RESOURCE_BARRIER::Transition(
        m_pipelineStatsBuffer.Get(),
        D3D12_RESOURCE_STATE_COPY_DEST,
        D3D12_RESOURCE_STATE_COPY_SOURCE);
    
    auto timestampBarrier = CD3DX12_RESOURCE_BARRIER::Transition(
        m_timestampBuffer.Get(),
        D3D12_RESOURCE_STATE_COPY_DEST,
        D3D12_RESOURCE_STATE_COPY_SOURCE);
    
    // Copy query data to readback buffers
    cmdList.ResourceBarrier(1, &pipelineStatsBarrier);
    cmdList.CopyResource(m_pipelineStatsReadbackBuffer.Get(), m_pipelineStatsBuffer.Get());
    
    cmdList.ResourceBarrier(1, &timestampBarrier);
    cmdList.CopyResource(m_timestampReadbackBuffer.Get(), m_timestampBuffer.Get());
       
    m_isCollecting = false;
}

void GPUStats::downloadStats()
{
    // Map the pipeline statistics readback buffer
    D3D12_QUERY_DATA_PIPELINE_STATISTICS* pPipelineStatsData = nullptr;
    CD3DX12_RANGE pipelineStatsRange(0, sizeof(D3D12_QUERY_DATA_PIPELINE_STATISTICS));
    ThrowIfFailed(m_pipelineStatsReadbackBuffer->Map(0, &pipelineStatsRange, reinterpret_cast<void**>(&pPipelineStatsData)));
    
    // Copy the pipeline statistics data
    m_downloadedStats = *pPipelineStatsData;
    
    // Unmap the buffer
    m_pipelineStatsReadbackBuffer->Unmap(0, nullptr);
    
    // Map the timestamp readback buffer
    uint64_t* pTimestampData = nullptr;
    CD3DX12_RANGE timestampRange(0, sizeof(uint64_t) * 2); // Space for begin and end timestamps
    ThrowIfFailed(m_timestampReadbackBuffer->Map(0, &timestampRange, reinterpret_cast<void**>(&pTimestampData)));
    
    // Copy the timestamp data
    m_timestampBegin = pTimestampData[0];
    m_timestampEnd = pTimestampData[1];
    
    // Unmap the buffer
    m_timestampReadbackBuffer->Unmap(0, nullptr);

    m_fDownloadedTimeMs = static_cast<double>(m_timestampEnd - m_timestampBegin) / static_cast<double>(m_timestampFrequency) * 1000.0;
}
 