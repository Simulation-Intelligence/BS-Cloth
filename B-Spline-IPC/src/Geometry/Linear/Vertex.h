#pragma once
#include "bsfwd.h"

namespace BSIPC
{
    struct SourceInfo
    {
    public:
        Float coeff;            // Coefficient at index of B-spline control point

        // When used as index of B-spline nodes, the index is the global index subject to the B-spline surface order. 
        // When used as index of vertex of linear triangle, the index is the global index in vector<Vec3> `mesh.vertices`
        UInt index;             // Index of B-spline or linear control point. 

        SourceInfo& operator= (const SourceInfo& sourceInfo)
        {
            this->coeff = sourceInfo.coeff;
            this->index = sourceInfo.index;
            return *this;
        }

        SourceInfo(const Float coeff, const UInt index) : coeff(coeff), index(index) {  }
        SourceInfo(const SourceInfo& sourceInfo) : coeff(sourceInfo.coeff), index(sourceInfo.index) {  }
    };

    struct Vertex
    {
    private:
        std::vector<SourceInfo> sources;

    public:
        Vec3 pos;

        const std::vector<SourceInfo>& GetSources() const { return this->sources; }

        Vertex& operator= (const Vertex& vertex)
        {
            this->pos = vertex.pos;
            this->sources = vertex.GetSources();
            return *this;
        }

        Vertex(const Vec3& pos, const std::vector<SourceInfo>& sources) : pos(pos), sources(sources) {  }
        Vertex(const Vertex& vertex) : pos(vertex.pos), sources(vertex.sources) {  }
    };
}
