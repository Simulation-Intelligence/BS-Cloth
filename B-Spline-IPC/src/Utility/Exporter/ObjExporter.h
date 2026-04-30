#pragma once
#include "bsfwd.h"

#include "Geometry/Linear/TriangularMesh.h"

namespace BSIPC
{
    class ObjExporter
    {
    private:
        std::optional<UInt> meshInterval = std::nullopt;
        std::vector<const TriangularMesh*> meshes;
        std::vector<Vec3> ctrlPts;

        Bool sepPatch, ctrl;
        std::vector<String> meshDirName;

        UInt SerializeMesh(std::ofstream& file, const TriangularMesh* mesh, const UInt offset) const;
        UInt SerializeHalfMesh(
            std::ofstream& file,
            const TriangularMesh* mesh, const UInt offset,
            const UInt step, Bool first
        ) const;
        UInt SerializeCube(std::ofstream& file, const Vec3 pos, const UInt offset, const Float edgeLength = 0.01) const;

        inline void SerializeHeader(std::ofstream& file) const
        {
            file << "# Exported by B-Spline-IPC\n";
            file << "s 1\n";            // enable smooth shading
        }

        inline static String SURFACE_W_STR = "surface_w";
        inline static String SURFACE_B_STR = "surface_b";
        inline static String CTRLS_STR = "ctrls";

    public:
        ObjExporter(Bool sepPatch, Bool ctrl, UInt meshCnt);

        inline void SetMeshInterval(const UInt interval) { this->meshInterval = interval; }
        inline void IncludeMesh(const TriangularMesh& mesh) { this->meshes.push_back(&mesh); }
        inline void IncludeCtrlPts(const std::vector<Vec3>& ctrlPts) { this->ctrlPts = ctrlPts; }
        inline void Clear() 
        { 
            this->meshes.clear();
            this->meshInterval = std::nullopt;
            this->ctrlPts.clear();
        } 

        void Serialize(const String filePrefix, UInt surfIndex, String frameIndexStr, Float cubeSize) const;
    };
}
