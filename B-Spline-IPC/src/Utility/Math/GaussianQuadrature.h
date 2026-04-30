#pragma once
#include "bsfwd.h"

#include <map>

namespace BSIPC
{
    struct QuadNodesData
    {
        std::vector<Float> nodes;
        std::vector<Float> weights;
    };

    class GaussianQuadrature
    {
    public:
        static const QuadNodesData GetNodesOfOrder(UInt order);

    private:
        // Stores the computed nodes into the cache
        static const QuadNodesData ComputePositiveNodes(UInt order);

    // Precomputed
    private:
        // All computed nodes should be in the interval [-1, 1]
        static std::map<UInt, QuadNodesData> cache;
        inline static const Float tolerance = 1e-10;
    };

    // Transform [-1, 1] canonic sample position for Gauss-Legendre to [DomainMin, DomainMax]
    Float TransformSamplePosition(Float sample, Float domainMin, Float domainMax);

    QuadNodesData TransformSamplePosition(const QuadNodesData& data, Float domainMin, Float domainMax);

    // Requires that the nodes are sorted
    void MergeConsecutiveNodes(std::vector<Float>& nodes, std::vector<Float>& weights, UInt ind);
}