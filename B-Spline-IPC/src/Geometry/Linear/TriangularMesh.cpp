#include "bspch.h"
#include "TriangularMesh.h"

namespace BSIPC
{
    void TriangularMesh::Clear()
    {
        this->vertices.clear();
        this->triangleIndices.clear();
    }

    void TriangularMesh::FixateRestPos()
    {
        this->restVertices.clear();
        for (const Vertex& vertex : this->vertices)
            this->restVertices.push_back(vertex.pos);
    }
}