#include "fem_parameters.h"
#include "IPCdistanceFuncs.h"
//#include <utility>
#include "oneapi/tbb.h"
#include "collisionUtil.h"

namespace IPC
{
    void SpatialHash::locateVoxelAxisIndex(const Vector3d& pos, Eigen::Array<int, 1, 3>& voxelAxisIndex)
    {
        voxelAxisIndex = ((pos - leftBottomCorner) * one_div_voxelSize).array().floor().template cast<int>();
    }

    int SpatialHash::locateVoxelIndex(const Vector3d& pos)
    {
        Eigen::Array<int, 1, 3> voxelAxisIndex;
        locateVoxelAxisIndex(pos, voxelAxisIndex);
        return voxelAxisIndex2VoxelIndex(voxelAxisIndex.data());
    }

    int SpatialHash::voxelAxisIndex2VoxelIndex(const int voxelAxisIndex[3])
    {
        return voxelAxisIndex2VoxelIndex(voxelAxisIndex[0], voxelAxisIndex[1], voxelAxisIndex[2]);
    }

    int SpatialHash::voxelAxisIndex2VoxelIndex(int ix, int iy, int iz)
    {
        return ix + iy * voxelCount[0] + iz * voxelCount0x1;
    }

    void SpatialHash::build(const  mesh3D& mesh, double voxelSize)
    {
        //double xmin = DBL_MAX, ymin = DBL_MAX, zmin = DBL_MAX;
        //double xmax = DBL_MIN, ymax = DBL_MIN, zmax = DBL_MIN;

        leftBottomCorner = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, (int)mesh.surfVerts.size()), Vector3d(1e32, 1e32, 1e32), [&](const tbb::blocked_range<int>& rg, Vector3d sceneSize) {
                for (int i = rg.begin(); i != rg.end(); i++) {
                    Vector3d pos = mesh.vertices[mesh.surfVerts[i]];
                    sceneSize[0] = min(pos[0], sceneSize[0]);
                    sceneSize[1] = min(pos[1], sceneSize[1]);
                    sceneSize[2] = min(pos[2], sceneSize[2]);
                }
                return sceneSize;
            },
            [&](Vector3d left, Vector3d right) {
                Vector3d output;
                output[0] = min(left[0], right[0]);
                output[1] = min(left[1], right[1]);
                output[2] = min(left[2], right[2]);
                return output;
            }
        );

        rightTopCorner = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, (int)mesh.surfVerts.size()), Vector3d(-1e32, -1e32, -1e32), [&](const tbb::blocked_range<int>& rg, Vector3d sceneSize) {
                for (int i = rg.begin(); i != rg.end(); i++) {
                    Vector3d pos = mesh.vertices[mesh.surfVerts[i]];
                    sceneSize[0] = max(pos[0], sceneSize[0]);
                    sceneSize[1] = max(pos[1], sceneSize[1]);
                    sceneSize[2] = max(pos[2], sceneSize[2]);
                }
                return sceneSize;
            },
            [&](Vector3d left, Vector3d right) {
                Vector3d output;
                output[0] = max(left[0], right[0]);
                output[1] = max(left[1], right[1]);
                output[2] = max(left[2], right[2]);
                return output;
            }
        );
        //for (int i = 0;i < mesh.surface.size();i++) {
        //    for (int j = 0;j < 3;j++) {
        //        Vector3d pos = mesh.vertices[mesh.surface[i][j]];
        //        if (xmin > pos[0]) xmin = pos[0];
        //        if (ymin > pos[1]) ymin = pos[1];
        //        if (zmin > pos[2]) zmin = pos[2];
        //        if (xmax < pos[0]) xmax = pos[0];
        //        if (ymax < pos[1]) ymax = pos[1];
        //        if (zmax < pos[2]) zmax = pos[2];
        //    }
        //}
        //leftBottomCorner = Vector3d(xmin, ymin, zmin);
        //rightTopCorner = Vector3d(xmax, ymax, zmax);

        one_div_voxelSize = 1.0 / voxelSize;
        Eigen::Array<double, 1, 3> range = rightTopCorner - leftBottomCorner;
        voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
        if (voxelCount.minCoeff() <= 0) {
            // cast overflow due to huge search direction
            one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
            voxelCount.setOnes();
        }
        voxelCount0x1 = voxelCount[0] * voxelCount[1];

        surfEdgeStartInd = mesh.surfVerts.size();
        surfTriStartInd = surfEdgeStartInd + mesh.surfEdges.size();

        std::vector<Eigen::Array<int, 1, 3>> svVoxelAxisIndex(mesh.surfVerts.size());
        std::vector<int> vI2SVI(mesh.vertexNum);

        tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI)
            //for (int svI = 0; svI < mesh.surfVerts.size(); ++svI)
            {
                int vI = mesh.surfVerts[svI];
                locateVoxelAxisIndex(mesh.vertices[vI], svVoxelAxisIndex[svI]);
                vI2SVI[vI] = svI;
            }
        );


        voxel.clear();


        for (int svI = 0; svI < mesh.surfVerts.size(); ++svI) {
            voxel[locateVoxelIndex(mesh.vertices[mesh.surfVerts[svI]])].emplace_back(svI);
        }

        std::vector<std::vector<int>> voxelLoc_e(mesh.surfEdges.size());

        tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int seCount)
            //for (int seCount = 0; seCount < mesh.surfEdges.size(); ++seCount)
            {
                const auto& seI = mesh.surfEdges[seCount];

                const Eigen::Array<int, 1, 3>& voxelAxisIndex_first = svVoxelAxisIndex[vI2SVI[seI.first]];
                const Eigen::Array<int, 1, 3>& voxelAxisIndex_second = svVoxelAxisIndex[vI2SVI[seI.second]];
                Eigen::Array<int, 1, 3> mins = voxelAxisIndex_first.min(voxelAxisIndex_second);
                Eigen::Array<int, 1, 3> maxs = voxelAxisIndex_first.max(voxelAxisIndex_second);
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            voxelLoc_e[seCount].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }
        );

        std::vector<std::vector<int>> voxelLoc_sf(mesh.surface.size());

        tbb::parallel_for(0, (int)mesh.surface.size(), 1, [&](int sfI)
            //for (int sfI = 0; sfI < mesh.surface.size(); ++sfI)
            {
                const Eigen::Array<int, 1, 3>& voxelAxisIndex0 = svVoxelAxisIndex[vI2SVI[mesh.surface[sfI][0]]];
                const Eigen::Array<int, 1, 3>& voxelAxisIndex1 = svVoxelAxisIndex[vI2SVI[mesh.surface[sfI][1]]];
                const Eigen::Array<int, 1, 3>& voxelAxisIndex2 = svVoxelAxisIndex[vI2SVI[mesh.surface[sfI][2]]];
                Eigen::Array<int, 1, 3> mins = voxelAxisIndex0.min(voxelAxisIndex1).min(voxelAxisIndex2);
                Eigen::Array<int, 1, 3> maxs = voxelAxisIndex0.max(voxelAxisIndex1).max(voxelAxisIndex2);
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            voxelLoc_sf[sfI].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }
        );

        for (int seCount = 0; seCount < voxelLoc_e.size(); ++seCount) {
            for (const auto& voxelI : voxelLoc_e[seCount]) {
                voxel[voxelI].emplace_back(seCount + surfEdgeStartInd);
            }
        }
        for (int sfI = 0; sfI < voxelLoc_sf.size(); ++sfI) {
            for (const auto& voxelI : voxelLoc_sf[sfI]) {
                voxel[voxelI].emplace_back(sfI + surfTriStartInd);
            }
        }
    }






    void SpatialHash::build(const  mesh3D& mesh, const vector<Vector3d>& searchDir, double& curMaxStepSize, double voxelSize, bool use_V_prev)
    {
        double pSize = 0;
        for (int svI = 0; svI < mesh.surfVerts.size(); ++svI) {
            int vI = mesh.surfVerts[svI];
            pSize += std::abs(searchDir[vI][0]);
            pSize += std::abs(searchDir[vI][1]);
            pSize += std::abs(searchDir[vI][2]);
        }
        pSize /= mesh.surfVerts.size() * 3;

        const double spanSize = curMaxStepSize * pSize / voxelSize;
        //if (spanSize > 1) {
        //    curMaxStepSize /= spanSize;
        //    // curMaxStepSize reduced for CCD spatial hash efficiency
        //}

        const vector<Vector3d>& V = use_V_prev ? mesh.V_prev : mesh.vertices;
        vector<Vector3d> SVt(mesh.surfVerts.size());
        std::vector<int> vI2SVI(mesh.vertexNum);

        tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI)
            //for (int svI = 0; svI < mesh.surfVerts.size(); ++svI) 
            {
                int vI = mesh.surfVerts[svI];
                SVt[svI] = V[vI] - curMaxStepSize * searchDir[vI];
                vI2SVI[vI] = svI;
            }
        );


        leftBottomCorner = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, (int)mesh.surfVerts.size()), Vector3d(1e32, 1e32, 1e32), [&](const tbb::blocked_range<int>& rg, Vector3d sceneSize) {
                for (int i = rg.begin(); i != rg.end(); i++) {
                    Vector3d pos = V[mesh.surfVerts[i]];
                    sceneSize[0] = min(pos[0], sceneSize[0]);
                    sceneSize[1] = min(pos[1], sceneSize[1]);
                    sceneSize[2] = min(pos[2], sceneSize[2]);

                    pos = SVt[i];
                    sceneSize[0] = min(pos[0], sceneSize[0]);
                    sceneSize[1] = min(pos[1], sceneSize[1]);
                    sceneSize[2] = min(pos[2], sceneSize[2]);

                }
                return sceneSize;
            },
            [&](Vector3d left, Vector3d right) {
                Vector3d output;
                output[0] = min(left[0], right[0]);
                output[1] = min(left[1], right[1]);
                output[2] = min(left[2], right[2]);
                return output;
            }
        );

        rightTopCorner = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, (int)mesh.surfVerts.size()), Vector3d(-1e32, -1e32, -1e32), [&](const tbb::blocked_range<int>& rg, Vector3d sceneSize) {
                for (int i = rg.begin(); i != rg.end(); i++) {
                    Vector3d pos = V[mesh.surfVerts[i]];
                    sceneSize[0] = max(pos[0], sceneSize[0]);
                    sceneSize[1] = max(pos[1], sceneSize[1]);
                    sceneSize[2] = max(pos[2], sceneSize[2]);

                    pos = SVt[i];
                    sceneSize[0] = max(pos[0], sceneSize[0]);
                    sceneSize[1] = max(pos[1], sceneSize[1]);
                    sceneSize[2] = max(pos[2], sceneSize[2]);
                }
                return sceneSize;
            },
            [&](Vector3d left, Vector3d right) {
                Vector3d output;
                output[0] = max(left[0], right[0]);
                output[1] = max(left[1], right[1]);
                output[2] = max(left[2], right[2]);
                return output;
            }
        );

        //double xmin = std::numeric_limits<double>::max(), ymin = std::numeric_limits<double>::max(), zmin = std::numeric_limits<double>::max();
        //double xmax = DBL_MIN, ymax = DBL_MIN, zmax = DBL_MIN;
        //for (int i = 0;i < mesh.surface.size();i++) {
        //    for (int j = 0;j < 3;j++) {
        //        Vector3d pos = V[mesh.surface[i][j]];
        //        if (xmin > pos[0]) xmin = pos[0];
        //        if (ymin > pos[1]) ymin = pos[1];
        //        if (zmin > pos[2]) zmin = pos[2];
        //        if (xmax < pos[0]) xmax = pos[0];
        //        if (ymax < pos[1]) ymax = pos[1];
        //        if (zmax < pos[2]) zmax = pos[2];

        //        pos = SVt[vI2SVI[mesh.surface[i][j]]];
        //        if (xmin > pos[0]) xmin = pos[0];
        //        if (ymin > pos[1]) ymin = pos[1];
        //        if (zmin > pos[2]) zmin = pos[2];
        //        if (xmax < pos[0]) xmax = pos[0];
        //        if (ymax < pos[1]) ymax = pos[1];
        //        if (zmax < pos[2]) zmax = pos[2];
        //    }
        //}

        //leftBottomCorner = Vector3d(xmin, ymin, zmin);
        //rightTopCorner = Vector3d(xmax, ymax, zmax);

        one_div_voxelSize = 1.0 / voxelSize;
        Eigen::Array<double, 1, 3> range = rightTopCorner - leftBottomCorner;
        voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
        if (voxelCount.minCoeff() <= 0) {
            // cast overflow due to huge search direction
            one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
            voxelCount.setOnes();
        }
        voxelCount0x1 = voxelCount[0] * voxelCount[1];

        surfEdgeStartInd = mesh.surfVerts.size();
        surfTriStartInd = surfEdgeStartInd + mesh.surfEdges.size();

        // precompute svVAI
        std::vector<Eigen::Array<int, 1, 3>> svMinVAI(mesh.surfVerts.size());
        std::vector<Eigen::Array<int, 1, 3>> svMaxVAI(mesh.surfVerts.size());

        tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI)
            {
                int vI = mesh.surfVerts[svI];
                Eigen::Array<int, 1, 3> v0VAI, vtVAI;
                locateVoxelAxisIndex(V[vI], v0VAI);
                locateVoxelAxisIndex(SVt[svI], vtVAI);
                svMinVAI[svI] = v0VAI.min(vtVAI);
                svMaxVAI[svI] = v0VAI.max(vtVAI);
            }

        );

        voxel.clear(); //TODO: parallel insert
        pointAndEdgeOccupancy.resize(0);
        pointAndEdgeOccupancy.resize(surfTriStartInd);


        tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI)

            {
                const Eigen::Array<int, 1, 3>& mins = svMinVAI[svI];
                const Eigen::Array<int, 1, 3>& maxs = svMaxVAI[svI];
                pointAndEdgeOccupancy[svI].reserve((maxs - mins + 1).prod());
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            pointAndEdgeOccupancy[svI].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }

        );



        tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int seCount)

            {
                int seIInd = seCount + surfEdgeStartInd;
                const auto& seI = mesh.surfEdges[seCount];

                Eigen::Array<int, 1, 3> mins = svMinVAI[vI2SVI[seI.first]].min(svMinVAI[vI2SVI[seI.second]]);
                Eigen::Array<int, 1, 3> maxs = svMaxVAI[vI2SVI[seI.first]].max(svMaxVAI[vI2SVI[seI.second]]);
                pointAndEdgeOccupancy[seIInd].reserve((maxs - mins + 1).prod());
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            pointAndEdgeOccupancy[seIInd].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }

        );

        std::vector<std::vector<int>> voxelLoc_sf(mesh.surface.size());

        tbb::parallel_for(0, (int)mesh.surface.size(), 1, [&](int sfI)

            {
                Eigen::Array<int, 1, 3> mins = svMinVAI[vI2SVI[mesh.surface[sfI][0]]].min(svMinVAI[vI2SVI[mesh.surface[sfI][1]]]).min(svMinVAI[vI2SVI[mesh.surface[sfI][2]]]);
                Eigen::Array<int, 1, 3> maxs = svMaxVAI[vI2SVI[mesh.surface[sfI][0]]].max(svMaxVAI[vI2SVI[mesh.surface[sfI][1]]]).max(svMaxVAI[vI2SVI[mesh.surface[sfI][2]]]);
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            voxelLoc_sf[sfI].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }

        );

        for (int i = 0; i < pointAndEdgeOccupancy.size(); ++i) {
            for (const auto& voxelI : pointAndEdgeOccupancy[i]) {
                voxel[voxelI].emplace_back(i);
            }
        }
        for (int sfI = 0; sfI < voxelLoc_sf.size(); ++sfI) {
            for (const auto& voxelI : voxelLoc_sf[sfI]) {
                voxel[voxelI].emplace_back(sfI + surfTriStartInd);
            }
        }
    }




    void SpatialHash::queryPointForTriangles(const Vector3d& pos, double radius, std::unordered_set<int>& triInds)
    {
        Eigen::Array<int, 1, 3> mins, maxs;
        locateVoxelAxisIndex(pos.array() - radius, mins);
        locateVoxelAxisIndex(pos.array() + radius, maxs);
        mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
        maxs = maxs.min(voxelCount.array() - 1);

        triInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfTriStartInd) {
                                triInds.insert(indI - surfTriStartInd);
                            }
                        }
                    }
                }
            }
        }
    }

    void SpatialHash::queryPointForTriangles(int svI, std::unordered_set<int>& sTriInds)
    {
        sTriInds.clear();
        for (const auto& voxelInd : pointAndEdgeOccupancy[svI]) {
            const auto& voxelI = voxel.find(voxelInd);
            assert(voxelI != voxel.end());
            for (const auto& indI : voxelI->second) {
                if (indI >= surfTriStartInd) {
                    sTriInds.insert(indI - surfTriStartInd);
                }
            }
        }
    }

    void SpatialHash::queryTriangleForPoints(const Vector3d& v0,
        const Vector3d& v1,
        const Vector3d& v2,
        double radius, std::unordered_set<int>& pointInds)
    {
        Vector3d leftBottom = v0.array().min(v1.array()).min(v2.array());
        Vector3d rightTop = v0.array().max(v1.array()).max(v2.array());
        Eigen::Array<int, 1, 3> mins, maxs;
        locateVoxelAxisIndex(leftBottom.array() - radius, mins);
        locateVoxelAxisIndex(rightTop.array() + radius, maxs);
        mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
        maxs = maxs.min(voxelCount - 1);

        pointInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI < surfEdgeStartInd) {
                                pointInds.insert(indI);
                            }
                        }
                    }
                }
            }
        }
    }



    void SpatialHash::queryTriangleForEdges(const Vector3d& v0,
        const Vector3d& v1,
        const Vector3d& v2,
        double radius, std::unordered_set<int>& edgeInds)
    {
        Vector3d leftBottom = v0.array().min(v1.array()).min(v2.array());
        Vector3d rightTop = v0.array().max(v1.array()).max(v2.array());
        Eigen::Array<int, 1, 3> mins, maxs;
        locateVoxelAxisIndex(leftBottom.array() - radius, mins);
        locateVoxelAxisIndex(rightTop.array() + radius, maxs);
        mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
        maxs = maxs.min(voxelCount - 1);

        edgeInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfEdgeStartInd && indI < surfTriStartInd) {
                                edgeInds.insert(indI - surfEdgeStartInd);
                            }
                        }
                    }
                }
            }
        }
    }

    void SpatialHash::queryEdgeForEdgesWithBBoxCheck(const mesh3D& mesh,
        const vector<Eigen::Vector3d>& searchDir, double curMaxStepSize,
        int seI, std::unordered_set<int>& sEdgeInds)
    {
        const Eigen::Matrix<double, 1, 3>& eI_v0 = mesh.vertices[mesh.surfEdges[seI].first];
        const Eigen::Matrix<double, 1, 3>& eI_v1 = mesh.vertices[mesh.surfEdges[seI].second];
        Eigen::Matrix<double, 1, 3> eI_v0t = eI_v0 - curMaxStepSize * searchDir[mesh.surfEdges[seI].first].transpose();
        Eigen::Matrix<double, 1, 3> eI_v1t = eI_v1 - curMaxStepSize * searchDir[mesh.surfEdges[seI].second].transpose();
        Eigen::Array<double, 1, 3> bboxEITopRight = eI_v0.array().max(eI_v0t.array()).max(eI_v1.array()).max(eI_v1t.array());
        Eigen::Array<double, 1, 3> bboxEIBottomLeft = eI_v0.array().min(eI_v0t.array()).min(eI_v1.array()).min(eI_v1t.array());
        sEdgeInds.clear();
        for (const auto& voxelInd : pointAndEdgeOccupancy[seI + surfEdgeStartInd]) {
            const auto& voxelI = voxel.find(voxelInd);
            assert(voxelI != voxel.end());
            for (const auto& indI : voxelI->second) {
                if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > seI) {
                    int seJ = indI - surfEdgeStartInd;
                    const Eigen::Matrix<double, 1, 3>& eJ_v0 = mesh.vertices[mesh.surfEdges[seJ].first];
                    const Eigen::Matrix<double, 1, 3>& eJ_v1 = mesh.vertices[mesh.surfEdges[seJ].second];
                    Eigen::Matrix<double, 1, 3> eJ_v0t = eJ_v0 - curMaxStepSize * searchDir[mesh.surfEdges[seJ].first].transpose();
                    Eigen::Matrix<double, 1, 3> eJ_v1t = eJ_v1 - curMaxStepSize * searchDir[mesh.surfEdges[seJ].second].transpose();
                    Eigen::Array<double, 1, 3> bboxEJTopRight = eJ_v0.array().max(eJ_v0t.array()).max(eJ_v1.array()).max(eJ_v1t.array());
                    Eigen::Array<double, 1, 3> bboxEJBottomLeft = eJ_v0.array().min(eJ_v0t.array()).min(eJ_v1.array()).min(eJ_v1t.array());
                    if (!((bboxEJBottomLeft - bboxEITopRight > 0.0).any() || (bboxEIBottomLeft - bboxEJTopRight > 0.0).any())) {
                        sEdgeInds.insert(indI - surfEdgeStartInd);
                    }
                }
            }
        }
    }


    void SpatialHash::queryEdgeForEdgesWithBBoxCheck(
        const mesh3D& mesh,
        const Eigen::Matrix<double, 1, 3>& vBegin,
        const Eigen::Matrix<double, 1, 3>& vEnd,
        double radius, std::vector<int>& edgeInds,
        int eIq)
    {
        // timer_mt.start(19);
        Eigen::Matrix<double, 1, 3> leftBottom = vBegin.array().min(vEnd.array()) - radius;
        Eigen::Matrix<double, 1, 3> rightTop = vBegin.array().max(vEnd.array()) + radius;
        Eigen::Array<int, 1, 3> mins, maxs;
        locateVoxelAxisIndex(leftBottom, mins);
        locateVoxelAxisIndex(rightTop, maxs);
        mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
        maxs = maxs.min(voxelCount.array() - 1);

        // timer_mt.start(20);
        edgeInds.resize(0);
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    // timer_mt.start(21);
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > eIq) {
                                int seJ = indI - surfEdgeStartInd;
                                const Eigen::Matrix<double, 1, 3>& eJ_v0 = mesh.vertices[mesh.surfEdges[seJ].first];
                                const Eigen::Matrix<double, 1, 3>& eJ_v1 = mesh.vertices[mesh.surfEdges[seJ].second];
                                Eigen::Array<double, 1, 3> bboxEJTopRight = eJ_v0.array().max(eJ_v1.array());
                                Eigen::Array<double, 1, 3> bboxEJBottomLeft = eJ_v0.array().min(eJ_v1.array());
                                if (!((bboxEJBottomLeft - rightTop.array() > 0.0).any() || (leftBottom.array() - bboxEJTopRight > 0.0).any())) {
                                    edgeInds.emplace_back(seJ);
                                }
                            }
                        }
                    }
                    // timer_mt.start(20);
                }
            }
        }
        std::sort(edgeInds.begin(), edgeInds.end());
        edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
        // timer_mt.stop();
    }

    void SpatialHash::queryEdgeForPE(const Vector3d& vBegin, const Vector3d& vEnd, std::vector<int>& svInds, std::vector<int>& edgeInds)
    {
        Vector3d leftBottom = vBegin.array().min(vEnd.array());
        Vector3d rightTop = vBegin.array().max(vEnd.array());
        Eigen::Array<int, 1, 3> mins, maxs;
        locateVoxelAxisIndex(leftBottom, mins);
        locateVoxelAxisIndex(rightTop, maxs);
        mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
        maxs = maxs.min(voxelCount - 1);

        // timer_mt.start(20);
        svInds.resize(0);
        edgeInds.resize(0);
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    // timer_mt.start(21);
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI < surfEdgeStartInd) {
                                svInds.emplace_back(indI);
                            }
                            else if (indI < surfTriStartInd) {
                                edgeInds.emplace_back(indI - surfEdgeStartInd);
                            }
                        }
                    }
                    // timer_mt.start(20);
                }
            }
        }
        std::sort(edgeInds.begin(), edgeInds.end());
        edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
        std::sort(svInds.begin(), svInds.end());
        svInds.erase(std::unique(svInds.begin(), svInds.end()), svInds.end());
    }

    void SpatialHash::queryPointForPrimitives(const Vector3d& pos, const Vector3d& dir, std::unordered_set<int>& sVInds, std::unordered_set<int>& sEdgeInds, std::unordered_set<int>& sTriInds)
    {
        Eigen::Array<int, 1, 3> mins, maxs;
        locateVoxelAxisIndex(pos.array().min((pos + dir).array()), mins);
        locateVoxelAxisIndex(pos.array().max((pos + dir).array()), maxs);
        mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
        maxs = maxs.min(voxelCount - 1);

        sVInds.clear();
        sEdgeInds.clear();
        sTriInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI < surfEdgeStartInd) {
                                sVInds.insert(indI);
                            }
                            else if (indI < surfTriStartInd) {
                                sEdgeInds.insert(indI - surfEdgeStartInd);
                            }
                            else {
                                sTriInds.insert(indI - surfTriStartInd);
                            }
                        }
                    }
                }
            }
        }
    }

    void SpatialHash::calculateActiveSet(mesh3D& mesh) {
        //mesh.Self_ActiveSet.resize(0);
        double sqrtDHat = std::sqrt(mesh.Hhat);
        vector<vector<MMCVID>> constraintSetPT(mesh.surfVerts.size());
        vector<vector<int>> cs_PT;
        cs_PT.resize(mesh.surfVerts.size());

        tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI)
            //for (int svI = 0; svI < mesh.surfVerts.size(); ++svI)
            {
                int vI = mesh.surfVerts[svI];
                std::unordered_set<int> triInds; //NOTE: different constraint order will result in numerically different results
                queryPointForTriangles(mesh.vertices[vI], sqrtDHat, triInds);
                for (const auto& sfI : triInds)
                {
                    const RowVector4i& sfVInd = mesh.surface[sfI].transpose();
                    if (!(mesh.boundaryTypes[vI] >= 2 && mesh.boundaryTypes[sfVInd[0]] >= 2 && mesh.boundaryTypes[sfVInd[1]] >= 2 && mesh.boundaryTypes[sfVInd[2]] >= 2))
                        if (!mesh.filterPointTrig(svI, sfI)) {
                            int dtype = dType_PT(mesh.vertices[vI], mesh.vertices[sfVInd[0]], mesh.vertices[sfVInd[1]], mesh.vertices[sfVInd[2]]);
                            double d;
                            switch (dtype) {
                            case 0: {
                                d_PP(mesh.vertices[vI], mesh.vertices[sfVInd[0]], d);
                                if (d < mesh.Hhat) {
                                    constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[0], -1, -1);
                                }
                                break;
                            }

                            case 1: {
                                d_PP(mesh.vertices[vI], mesh.vertices[sfVInd[1]], d);
                                if (d < mesh.Hhat) {
                                    constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[1], -1, -1);
                                }
                                break;
                            }

                            case 2: {
                                d_PP(mesh.vertices[vI], mesh.vertices[sfVInd[2]], d);
                                if (d < mesh.Hhat) {
                                    constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[2], -1, -1);
                                }
                                break;
                            }

                            case 3: {
                                d_PE(mesh.vertices[vI], mesh.vertices[sfVInd[0]], mesh.vertices[sfVInd[1]], d);
                                if (d < mesh.Hhat) {
                                    constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[0], sfVInd[1], -1);
                                }
                                break;
                            }

                            case 4: {
                                d_PE(mesh.vertices[vI], mesh.vertices[sfVInd[1]], mesh.vertices[sfVInd[2]], d);
                                if (d < mesh.Hhat) {
                                    constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[1], sfVInd[2], -1);
                                }
                                break;
                            }

                            case 5: {
                                d_PE(mesh.vertices[vI], mesh.vertices[sfVInd[2]], mesh.vertices[sfVInd[0]], d);
                                if (d < mesh.Hhat) {
                                    constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[2], sfVInd[0], -1);
                                }
                                break;
                            }

                            case 6: {
                                d_PT(mesh.vertices[vI], mesh.vertices[sfVInd[0]], mesh.vertices[sfVInd[1]], mesh.vertices[sfVInd[2]], d);
                                if (d < mesh.Hhat) {
                                    constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[0], sfVInd[1], sfVInd[2]);
                                }
                                break;
                            }

                            default:
                                break;
                            }
                            //d_PE(mesh.vertices[2013], mesh.vertices[1930], mesh.vertices[632], d);
                            //if (d < 4.3e-6) {
                            //printf("%f\n", d);
                            //}
                            cs_PT[svI].emplace_back(sfI);
                        }
                }
            }
        );
        vector<vector<MMCVID>> constraintSetEE(mesh.surfEdges.size());
        vector<vector<int>> cs_EE;
        cs_EE.resize(mesh.surfEdges.size());

        tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int eI)
            //for (int eI = 0; eI < mesh.surfEdges.size(); ++eI)
            {
                const auto& meshEI = mesh.surfEdges[eI];
                vector<int> edgeInds; //NOTE: different constraint order will result in numerically different results
                queryEdgeForEdgesWithBBoxCheck(mesh, mesh.vertices[meshEI.first].transpose(), mesh.vertices[meshEI.second].transpose(), sqrtDHat, edgeInds, eI);
                for (const auto& eJ : edgeInds) {
                    const auto& meshEJ = mesh.surfEdges[eJ];
                    if (!(mesh.boundaryTypes[meshEI.first] >= 2 && mesh.boundaryTypes[meshEI.second] >= 2 && mesh.boundaryTypes[meshEJ.first] >= 2 && mesh.boundaryTypes[meshEJ.second] >= 2))
                        if (!(mesh.filterEdgeEdge(eI, eJ) || eI > eJ)) {

                            int dtype = dType_EE(mesh.vertices[meshEI.first], mesh.vertices[meshEI.second], mesh.vertices[meshEJ.first], mesh.vertices[meshEJ.second]);

                            double EECrossSqNorm, eps_x;
                            computeEECrossSqNorm(mesh.vertices[meshEI.first], mesh.vertices[meshEI.second], mesh.vertices[meshEJ.first], mesh.vertices[meshEJ.second], EECrossSqNorm);
                            compute_eps_x(mesh, meshEI.first, meshEI.second, meshEJ.first, meshEJ.second, eps_x);

                            int add_e = (EECrossSqNorm < eps_x) ? -eJ - 2 : -1;
                            double d;
                            switch (dtype) {
                            case 0: {
                                d_PP(mesh.vertices[meshEI.first], mesh.vertices[meshEJ.first], d);
                                if (d < mesh.Hhat) {
                                    constraintSetEE[eI].emplace_back(-meshEI.first - 1, meshEJ.first, -1, add_e);
                                }
                                break;
                            }

                            case 1: {
                                d_PP(mesh.vertices[meshEI.first], mesh.vertices[meshEJ.second], d);
                                if (d < mesh.Hhat) {
                                    constraintSetEE[eI].emplace_back(-meshEI.first - 1, meshEJ.second, -1, add_e);
                                }
                                break;
                            }

                            case 2: {
                                d_PE(mesh.vertices[meshEI.first], mesh.vertices[meshEJ.first], mesh.vertices[meshEJ.second], d);
                                if (d < mesh.Hhat) {
                                    constraintSetEE[eI].emplace_back(-meshEI.first - 1, meshEJ.first, meshEJ.second, add_e);
                                }
                                break;
                            }

                            case 3: {
                                d_PP(mesh.vertices[meshEI.second], mesh.vertices[meshEJ.first], d);
                                if (d < mesh.Hhat) {
                                    constraintSetEE[eI].emplace_back(-meshEI.second - 1, meshEJ.first, -1, add_e);
                                }
                                break;
                            }

                            case 4: {
                                d_PP(mesh.vertices[meshEI.second], mesh.vertices[meshEJ.second], d);
                                if (d < mesh.Hhat) {
                                    constraintSetEE[eI].emplace_back(-meshEI.second - 1, meshEJ.second, -1, add_e);
                                }
                                break;
                            }

                            case 5: {
                                d_PE(mesh.vertices[meshEI.second], mesh.vertices[meshEJ.first], mesh.vertices[meshEJ.second], d);
                                if (d < mesh.Hhat) {
                                    constraintSetEE[eI].emplace_back(-meshEI.second - 1, meshEJ.first, meshEJ.second, add_e);
                                }
                                break;
                            }

                            case 6: {
                                d_PE(mesh.vertices[meshEJ.first], mesh.vertices[meshEI.first], mesh.vertices[meshEI.second], d);
                                if (d < mesh.Hhat) {
                                    constraintSetEE[eI].emplace_back(-meshEJ.first - 1, meshEI.first, meshEI.second, add_e);
                                }
                                break;
                            }

                            case 7: {
                                d_PE(mesh.vertices[meshEJ.second], mesh.vertices[meshEI.first], mesh.vertices[meshEI.second], d);
                                if (d < mesh.Hhat) {
                                    constraintSetEE[eI].emplace_back(-meshEJ.second - 1, meshEI.first, meshEI.second, add_e);
                                }
                                break;
                            }

                            case 8: {
                                d_EE(mesh.vertices[meshEI.first], mesh.vertices[meshEI.second], mesh.vertices[meshEJ.first], mesh.vertices[meshEJ.second], d);
                                if (d < mesh.Hhat) {
                                    if (add_e <= -2) {
                                        constraintSetEE[eI].emplace_back(meshEI.first, meshEI.second, meshEJ.first, -meshEJ.second - mesh.surfEdges.size() - 2);
                                    }
                                    else {
                                        constraintSetEE[eI].emplace_back(meshEI.first, meshEI.second, meshEJ.first, meshEJ.second);
                                    }
                                }
                                break;
                            }

                            default:
                                break;
                            }


                            cs_EE[eI].emplace_back(eJ);

                        }

                }
            }
        );

        mesh.Self_CCD_ActiveSet.resize(0);
        mesh.Self_CCD_ActiveSet.reserve(cs_PT.size() + cs_EE.size());
        for (int svI = 0; svI < cs_PT.size(); ++svI) {
            for (const auto& sfI : cs_PT[svI]) {
                mesh.Self_CCD_ActiveSet.emplace_back(-svI - 1, sfI);
            }
        }
        cs_PT.clear();
        for (int eI = 0; eI < cs_EE.size(); ++eI) {
            for (const auto& eJ : cs_EE[eI]) {
                mesh.Self_CCD_ActiveSet.emplace_back(eI, eJ);
            }
        }
        cs_EE.clear();

        mesh.Self_ActiveSet.resize(0);
        mesh.Self_ActiveSet.reserve(constraintSetPT.size() + constraintSetEE.size());

        std::map<MMCVID, int> constraintCounter;
        for (const auto& csI : constraintSetPT) {
            for (const auto& cI : csI) {
                if (cI.data[3] < 0) {
                    // PP or PE
                    ++constraintCounter[cI];
                }
                else {
                    mesh.Self_ActiveSet.emplace_back(cI);
                }
            }
        }
        constraintSetPT.clear();
        mesh.Self_EE_ActiveSet.resize(0);
        mesh.Self_EEeIe_ActiveSet.resize(0);
        int eI = 0;
        for (const auto& csI : constraintSetEE) {
            for (const auto& cI : csI) {
                if (cI.data[3] >= 0) {
                    // regular EE
                    mesh.Self_ActiveSet.emplace_back(cI);
                }
                else if (cI.data[3] == -1) {
                    // regular PP or PE
                    ++constraintCounter[cI];
                }
                else if (cI.data[3] >= -mesh.surfEdges.size() - 1) {
                    // nearly parallel PP or PE
                    mesh.Self_EE_ActiveSet.emplace_back(cI.data[0], cI.data[1], cI.data[2], -1);
                    mesh.Self_EEeIe_ActiveSet.emplace_back(eI, -cI.data[3] - 2);
                }
                else {
                    // nearly parallel EE
                    mesh.Self_EE_ActiveSet.emplace_back(cI.data[0], cI.data[1], cI.data[2], -cI.data[3] - mesh.surfEdges.size() - 2);
                    mesh.Self_EEeIe_ActiveSet.emplace_back(-1, -1);
                }
            }
            ++eI;
        }
        constraintSetEE.clear();

        mesh.Self_ActiveSet.reserve(mesh.Self_ActiveSet.size() + constraintCounter.size());
        for (const auto& ccI : constraintCounter) {
            mesh.Self_ActiveSet.emplace_back(ccI.first.data[0], ccI.first.data[1], ccI.first.data[2], -ccI.second);
        }
    }


    void Ground::Init(const Vector3d& m_normal, const double& d) {
        normal = m_normal;
        D = d;
    }

    Ground::Ground() {
        Init(Vector3d(0, 1, 0), -1);
    }

    double Ground::calculateGapFromObj(const mesh3D& mesh, const int& vId) const {
        double dist = normal.dot(mesh.vertices[vId]) - D;
        return dist * dist;
    }
    
    //double Ground::calculateConstraintVal(const mesh3D& mesh, const int& vId) {
    //    double dist = normal.dot(mesh.vertices[vId]) - D;
    //    return dist * dist;
    //}

    void Ground::calculateActiveSet(mesh3D& mesh) const {
        mesh.Environment_ActiveSet.resize(0);


        tbb::spin_mutex groundMutex;//, countMutex3;
        tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int i)
            //for (int i = 0;i < mesh.surfVerts.size();i++) 
            {
                double dis = calculateGapFromObj(mesh, mesh.surfVerts[i]);
                if (dis < mesh.Hhat) {
                    groundMutex.lock();
                    mesh.Environment_ActiveSet.push_back(mesh.surfVerts[i]);
                    groundMutex.unlock();
                }
            }
        );
    }
}