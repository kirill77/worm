#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <d3d12.h>
#include <wrl/client.h>

// Forward declarations
namespace Microsoft { namespace WRL { template<typename> class ComPtr; } }
struct ID3D12Device;
struct ID3D12QueryHeap;
struct ID3D12Resource;
class GPUQueue;

class GPUStats
{
public:
    GPUStats(Microsoft::WRL::ComPtr<ID3D12Device> device, std::shared_ptr<GPUQueue> gpuQueue);
    ~GPUStats();

    // Begin collecting GPU statistics
    void begin();

    // End collecting GPU statistics
    void end();

    // Get the collected statistics as a string
    std::string getStats() const;

private:
    // Initialize query resources
    void initializeQueries();

    // Read back query data
    void readQueryData();

    Microsoft::WRL::ComPtr<ID3D12Device> m_device;
    std::shared_ptr<GPUQueue> m_gpuQueue;
    Microsoft::WRL::ComPtr<ID3D12QueryHeap> m_queryHeap;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_queryBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_queryReadbackBuffer;
    
    // Query types we're tracking
    enum class QueryType
    {
        PipelineStatistics,
        Timestamp
    };
    
    // Query indices
    struct QueryIndices
    {
        uint32_t pipelineStatsBegin;
        uint32_t pipelineStatsEnd;
        uint32_t timestampBegin;
        uint32_t timestampEnd;
    } m_queryIndices;
    
    // Pipeline statistics data
    D3D12_QUERY_DATA_PIPELINE_STATISTICS m_pipelineStats;
    
    // Timestamp data
    uint64_t m_timestampBegin;
    uint64_t m_timestampEnd;
    
    // Timestamp frequency for converting to milliseconds
    uint64_t m_timestampFrequency;
    
    // Flag to indicate if we're currently collecting stats
    bool m_isCollecting;
}; 