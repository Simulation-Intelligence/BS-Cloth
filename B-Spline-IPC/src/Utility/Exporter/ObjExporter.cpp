#include "bspch.h"
#include "ObjExporter.h"

#include <filesystem>

#include "Utility/OS/Command.h"

namespace BSIPC
{
    namespace fs = std::filesystem;

    ObjExporter::ObjExporter(Bool sepPatch, Bool ctrl, UInt meshCnt) : sepPatch(sepPatch), ctrl(ctrl)
    {
        fs::path objPath = fs::path("objs");

        this->meshDirName.resize(meshCnt);

        for (UInt i = 0; i != meshCnt; ++i)
        {
            String surfStr = "surf_" + std::to_string(i);
            fs::path meshPath = objPath / surfStr;
            OS::CreateDir(meshPath.string());
            this->meshDirName[i] = surfStr;

            if (sepPatch)
            {
                fs::path surfaceWPath = meshPath / ObjExporter::SURFACE_W_STR;
                fs::path surfaceBPath = meshPath / ObjExporter::SURFACE_B_STR;
                OS::CreateDir(surfaceWPath.string());
                OS::CreateDir(surfaceBPath.string());
            }

            if (ctrl)
            {
                fs::path ctrlsPath = meshPath / ObjExporter::CTRLS_STR;
                OS::CreateDir(ctrlsPath.string());
            }
        }
    }

    UInt ObjExporter::SerializeMesh(std::ofstream& file, TriangularMesh const* mesh, const UInt offset) const
    {
        const std::vector<Vertex>& vertices = mesh->GetVertices();
        const std::vector<Vec3>& normals = mesh->GetNormals();
        const std::vector<IVec3>& triangleIndices = mesh->GetTriangleIndices();

        // List all the vertices
        for (size_t i = 0; i != vertices.size(); ++i)
            file << "v " << vertices[i].pos[0] << " " << vertices[i].pos[1] << " " << vertices[i].pos[2] << "\n";

        // List all the normals
        for (size_t i = 0; i != normals.size(); ++i)
            file << "vn " << normals[i][0] << " " << normals[i][1] << " " << normals[i][2] << "\n";

        // Compose all faces. Obj face indices starts with 1
        for (size_t i = 0; i != triangleIndices.size(); ++i)
        {
            std::array<UInt, 3> indices = {
                triangleIndices[i][0] + offset,
                triangleIndices[i][1] + offset,
                triangleIndices[i][2] + offset
            };

            String faceStr = fmt::format("f {}//{} {}//{} {}//{}\n",
                indices[0], indices[0],
                indices[1], indices[1],
                indices[2], indices[2]
            );

            file << faceStr;
        }
        
        return offset + static_cast<UInt>(vertices.size());
    }

    UInt ObjExporter::SerializeHalfMesh(
        std::ofstream& file, 
        const TriangularMesh* mesh, const UInt offset, 
        const UInt step, Bool first
    ) const
    {
        const std::vector<Vertex>& vertices = mesh->GetVertices();
        const std::vector<IVec3>& triangleIndices = mesh->GetTriangleIndices();

        // List all the vertices
        for (size_t i = 0; i != vertices.size(); ++i)
            file << "v " << vertices[i].pos[0] << " " << vertices[i].pos[1] << " " << vertices[i].pos[2] << "\n";

        // Compose all faces. Obj face indices starts with 1
        UInt facesAdded = 0;
        UInt modulo = first ? 0 : 1;
        Bool evenSideGrid = mesh->GetTriangleIndices().size() / step % 2 == 0;
        UInt sideLength = static_cast<UInt>(std::sqrt(mesh->GetTriangleIndices().size() / step));
        for (size_t i = 0; i != triangleIndices.size(); ++i)
            if ((!evenSideGrid && (i / step) % 2 == modulo) || (evenSideGrid && (i / step + i / (step * sideLength)) % 2 == modulo))
            {
                file << "f " << triangleIndices[i][0] + offset << " " << triangleIndices[i][1] + offset << " " << triangleIndices[i][2] + offset << "\n";
                ++facesAdded;
            }

        return offset + facesAdded;
    }

    void ObjExporter::Serialize(const String filePrefix, UInt surfIndex, String frameIndexStr, Float cubeSize = 0.01) const
    {
        BSIPC_ASSERT(this->meshInterval.has_value() == this->sepPatch, "patch separation info should be consistent");
        BSIPC_ASSERT(this->ctrlPts.empty() == !this->ctrl, "ctrlPts info should be consistent");

        fs::path surfPath = fs::path("objs") / this->meshDirName[surfIndex];

        if (this->meshInterval == std::nullopt)
        {
            String filename = fmt::format("{}{}.obj", filePrefix, frameIndexStr);
            fs::path meshPath = surfPath / filename;
            std::ofstream file(meshPath.string());
            this->SerializeHeader(file);

            UInt curFaceOffset = 1;

            // Serialize all the meshes
            for (UInt i = 0; i != this->meshes.size(); ++i)
            {
                file << "o mesh" << i << "\n";
                curFaceOffset = this->SerializeMesh(file, this->meshes[i], curFaceOffset);
            }

            file.close();
        }
        else
        {
            std::array<String, 2> dirs = { ObjExporter::SURFACE_W_STR, ObjExporter::SURFACE_B_STR };
            std::array<String, 2> filenames = 
            { 
                fmt::format("{}w_{}.obj", filePrefix, frameIndexStr),
                fmt::format("{}b_{}.obj", filePrefix, frameIndexStr)
            };
            std::array<fs::path, 2> paths = { surfPath / dirs[0] / filenames[0], surfPath / dirs[1] / filenames[1] };
            std::array<std::ofstream, 2> files = { std::ofstream(paths[0].string()), std::ofstream(paths[1].string())};
            std::array<UInt, 2> offsets = { 1, 1 };
            for (auto& file : files)
            {
                this->SerializeHeader(file);
                file << "o mesh" << "\n";
            }

            for (UInt i = 0; i != this->meshes.size(); ++i)
            {
                offsets[0] = this->SerializeHalfMesh(files[0], this->meshes[i], offsets[0], this->meshInterval.value(), true);
                offsets[1] = this->SerializeHalfMesh(files[1], this->meshes[i], offsets[1], this->meshInterval.value(), false);
            }

            for (auto& file : files)
                file.close();
        }

        if (!this->ctrlPts.empty())
        {
            String filename = fmt::format("{}ctrlPts_{}.obj", filePrefix, frameIndexStr);
            fs::path meshPath = surfPath / fs::path(ObjExporter::CTRLS_STR) / filename;
            std::ofstream file(meshPath.string());
            this->SerializeHeader(file);
            UInt curFaceOffset = 1;
            for (UInt i = 0; i != this->ctrlPts.size(); ++i)
            {
                file << "o ctrlPt" << i << "\n";
                curFaceOffset = this->SerializeCube(file, this->ctrlPts[i], curFaceOffset, cubeSize);
            }
            file.close();
        }
    }

    UInt ObjExporter::SerializeCube(std::ofstream& file, const Vec3 pos, const UInt offset, const Float edgeLength) const
    {
        const Float halfEdge = edgeLength / 2.0f;
        const Vec3 vertices[8] = {
            Vec3(pos[0] - halfEdge, pos[1] - halfEdge, pos[2] - halfEdge),
            Vec3(pos[0] + halfEdge, pos[1] - halfEdge, pos[2] - halfEdge),
            Vec3(pos[0] + halfEdge, pos[1] + halfEdge, pos[2] - halfEdge),
            Vec3(pos[0] - halfEdge, pos[1] + halfEdge, pos[2] - halfEdge),
            Vec3(pos[0] - halfEdge, pos[1] - halfEdge, pos[2] + halfEdge),
            Vec3(pos[0] + halfEdge, pos[1] - halfEdge, pos[2] + halfEdge),
            Vec3(pos[0] + halfEdge, pos[1] + halfEdge, pos[2] + halfEdge),
            Vec3(pos[0] - halfEdge, pos[1] + halfEdge, pos[2] + halfEdge)
        };

        const IVec3 faces[12] = 
        {
            IVec3(0, 1, 2), IVec3(0, 2, 3),
            IVec3(4, 6, 5), IVec3(4, 7, 6),
            IVec3(0, 4, 1), IVec3(1, 4, 5),
            IVec3(1, 5, 2), IVec3(2, 5, 6),
            IVec3(2, 6, 3), IVec3(3, 6, 7),
            IVec3(3, 7, 0), IVec3(0, 7, 4)
        };

        for (UInt i = 0; i != 8; ++i)
            file << "v " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << "\n";

        for (UInt i = 0; i != 12; ++i)
            file << "f " << faces[i][0] + offset << " " << faces[i][1] + offset << " " << faces[i][2] + offset << "\n";

        return offset + 8;
    }
}
