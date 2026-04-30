#include "atrpch.h"
#include "Mesh.h"

namespace ATR
{
    void Mesh::AddTriangle(std::array<Vertex, 3> vertices)
    {
        for (auto& vertex : vertices)
        {
            auto foundIter = std::find(this->vertices.begin(), this->vertices.end(), vertex);
            UInt newIndex = 0;
            if (foundIter == this->vertices.cend())
            {
                UInt size = static_cast<UInt>(this->vertices.size());
                newIndex = size;
                this->vertices.push_back(vertex);
            }
            else
                newIndex = static_cast<UInt>(foundIter - this->vertices.begin());

            this->indices.push_back(newIndex);
        }
    }

    void Mesh::UpdateMesh(const Mesh& mesh)
    {
        this->Clear();

        this->indices = mesh.GetIndices();
        this->vertices = mesh.GetVertices();
    }

    String Mesh::Info() const
    {
        String msg;

        for (const auto& vertex : this->vertices)
            msg += vertex.Info() + "\n";

        return msg;
    }
}
