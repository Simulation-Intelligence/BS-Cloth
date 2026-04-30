#pragma once

#include "atrfwd.h"

#include "Vertex.h"

namespace ATR
{
    class Mesh
    {
    public:
        void AddTriangle(std::array<Vertex, 3> vertices);

        inline void FormTriangle(UInt x, UInt y, UInt z) { this->indices.push_back(x); this->indices.push_back(y); this->indices.push_back(z); }
        inline void SetVertices(std::vector<Vertex>& vertices) { this->vertices = vertices; }

        inline const std::vector<UInt>& GetIndices() const { return this->indices; }
        inline const std::vector<Vertex>& GetVertices() const { return this->vertices; }

        inline void UpdateVertex(UInt index, Vec3 pos, Vec3 normal) { this->vertices[index].pos = pos; this->vertices[index].normal = normal; }

        inline void Clear() { this->indices.clear(); this->vertices.clear(); }

        void UpdateMesh(const Mesh& mesh);

        String Info() const;

    private:
        std::vector<UInt> indices;
        std::vector<Vertex> vertices;
    };

}
