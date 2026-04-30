#pragma once
#ifndef _COLLISTION_UTIL_H_
#define _COLLISTION_UTIL_H_

#include "mesh.h"
#include <unordered_map>
#include <unordered_set>

namespace IPC
{

    using namespace Eigen;

    class SpatialHash {
    public: // data
        Vector3d leftBottomCorner, rightTopCorner;
        double one_div_voxelSize;
        Array<int, 1, 3> voxelCount;
        int voxelCount0x1;

        int surfEdgeStartInd, surfTriStartInd;

        std::unordered_map<int, std::vector<int>> voxel;
        std::vector<std::vector<int>> pointAndEdgeOccupancy;

    public: // constructor
        SpatialHash(void) {}
        SpatialHash(const mesh3D& mesh, double voxelSize, bool use_V_prev = false);
        SpatialHash(const mesh3D& mesh, const Eigen::VectorXd& searchDir, double curMaxStepSize, double voxelSize, bool use_V_prev = false);

        void build(const mesh3D& mesh, double voxelSize);
        void build(const mesh3D& mesh, const vector<Vector3d>& searchDir, double& curMaxStepSize, double voxelSize, bool use_V_prev = false);

        void locateVoxelAxisIndex(const Vector3d& pos, Eigen::Array<int, 1, 3>& voxelAxisIndex);
        int locateVoxelIndex(const Vector3d& pos);
        int voxelAxisIndex2VoxelIndex(const int voxelAxisIndex[3]);
        int voxelAxisIndex2VoxelIndex(int ix, int iy, int iz);
        // finds triangles intersecting the cube [[ pos - radius, pos + radius ]]
        void queryPointForTriangles(const Vector3d& pos, double radius, std::unordered_set<int>& triInds);
        void queryPointForTriangles(int svI, std::unordered_set<int>& sTriInds);
        void queryTriangleForPoints(const Vector3d& v0,
            const Vector3d& v1,
            const Vector3d& v2,
            double radius, std::unordered_set<int>& pointInds);
        void queryTriangleForEdges(const Vector3d& v0,
            const Vector3d& v1,
            const Vector3d& v2,
            double radius, std::unordered_set<int>& edgeInds);
        void queryPointForPrimitives(const Vector3d& pos,
            const Vector3d& dir, std::unordered_set<int>& sVInds,
            std::unordered_set<int>& sEdgeInds, std::unordered_set<int>& sTriInds);
        void queryEdgeForPE(const Vector3d& vBegin,
            const Vector3d& vEnd,
            std::vector<int>& svInds, std::vector<int>& edgeInds);
        void queryEdgeForEdgesWithBBoxCheck(
            const mesh3D& mesh,
            const Eigen::Matrix<double, 1, 3>& vBegin,
            const Eigen::Matrix<double, 1, 3>& vEnd,
            double radius, std::vector<int>& edgeInds,
            int eIq = -1);
        void queryEdgeForEdgesWithBBoxCheck(const mesh3D& mesh,
            const vector<Eigen::Vector3d>& searchDir, double curMaxStepSize,
            int seI, std::unordered_set<int>& sEdgeInds);
        void calculateActiveSet(mesh3D& mesh);
    };

    class Ground {
    public: // data
        Vector3d normal;
        double D;
    public:
        Ground();
        void Init(const Vector3d& m_normal, const double& d);
        double calculateGapFromObj(const mesh3D& mesh, const int& vId) const;
        void calculateActiveSet(mesh3D& mesh) const;
        //void calculateConstraintVal();
    };


}

#endif // ! _COLLISTION_UTIL_H_

