#pragma once
#include "bsfwd.h"

#include "Vertex.h"

namespace BSIPC
{
    class TriangularMesh
    {
    private:
        std::vector<Vertex> vertices;
        std::vector<Vec3> normals;
        std::vector<Vec3> restVertices;

        /// The following objects are all described using indices in the vertices vector
        std::vector<IVec2> edgeIndices;
        std::vector<IVec3> triangleIndices;

    public:
        TriangularMesh() = default;

        /* ***************************** Getter/Setter ***************************** */

        void Clear();
        void FixateRestPos();
         
        inline const std::vector<Vertex>& GetVertices() const { return this->vertices; }
        inline const std::vector<Vec3>& GetNormals() const { return this->normals; }
        inline const std::vector<Vec3>& GetRestVerticesPos() const { return this->restVertices; }

        inline const std::vector<IVec2>& GetEdgeIndices() const { return this->edgeIndices; }
        inline const std::vector<IVec3>& GetTriangleIndices() const { return this->triangleIndices; }

        inline void SetVertices(const std::vector<Vertex>& vertices) { this->vertices = vertices; }
        inline void SetNormals(const std::vector<Vec3>& normals) { this->normals = normals; }
        inline void UpdateVertex(size_t index, const Vec3 pos) { this->vertices[index].pos = pos; }
        inline void UpdateNormal(size_t index, const Vec3 normal) { this->normals[index] = normal; }
        inline void FormTriangle(const UInt x, const UInt y, const UInt z) 
        { 
            this->triangleIndices.emplace_back(x, y, z);
            this->edgeIndices.emplace_back(x, y);
            this->edgeIndices.emplace_back(y, z);
            this->edgeIndices.emplace_back(z, x);
        }
    };
}
