#include "bspch.h"

#include "Solver.h"

#include <fstream>

#include "Utility/Math/MathHelpers.h"
#include "Utility/Loader/YAMLLoader.h"
#include "Utility/EigenHelper.h"
#include "Utility/Profiler.h"

#include "IPC/IPC_FUNC.h"


namespace BSIPC
{
    Solver::Solver(SolverConfig config) :
        config(config),
        contactStiffness(0.),
        trigVertCntSum(0),
        bsVertCntSum(0),
        bsPatchCntSum(0),
        inertiaHessEntryCnt(0),
        potHessEntryCnt(0),
        curIndexInHessCache(0),
        simlTime(0),
        seamConstraintCnt(0)
    {  
        this->ground.Init(Vec3(0, 0, 1), config.groundHeight);
        if (config.contactStiffness.has_value())
            this->contactStiffness = config.contactStiffness.value();

#if not defined BSIPC_DISABLE_RENDERER
        this->renderer = nullptr;
#endif
    }

#if defined BSIPC_NUMERIC_TEST
    DMat Solver::NumericalBarrierGrad_Trig()
    {
        UInt size = this->trigVertCntSum;

        auto BarrierVal = [&]() -> Float
        {
            IPC::BuildCollisionSets(this->contactMesh, this->spatialHash, this->ground, true);
            Float val = IPC::ComputeBarrierEnergyVal(this->contactMesh, this->ground, this->contactStiffness);
            return val;
        };

        DMat grad = DMat::Zero(3 * size, 1);
        Float h = 1e-8;
        for (UInt i = 0; i != size; ++i)
        {
            for (UInt indexI = 0; indexI != 3; ++indexI)
            {
                Vec3& perturbVertex = this->contactMesh.vertices[i];
                perturbVertex[indexI] += h;
                Float upVal = BarrierVal();

                perturbVertex[indexI] -= 2 * h;
                Float downVal = BarrierVal();

                perturbVertex[indexI] += h;
                grad(3 * i + indexI, 0) = (upVal - downVal) / (2 * h);
            }
        }

        return grad;
    }

    DMat Solver::NumericalBarrierHess_Trig()
    {
        UInt size = this->trigVertCntSum;

        auto BarrierVal = [&]() -> Float
        {
            IPC::BuildCollisionSets(this->contactMesh, this->spatialHash, this->ground, true);
            Float val = IPC::ComputeBarrierEnergyVal(this->contactMesh, this->ground, this->contactStiffness);
            return val;
        };

        DMat hess = DMat::Zero(3 * size, 3 * size);
        Float h = 1e-8;
        Float hSq = h * h;

        for (UInt i = 0; i != size; ++i)
        {
            for (UInt j = 0; j != size; ++j)
            {
                for (UInt indexI = 0; indexI != 3; ++indexI)
                {
                    for (UInt indexJ = 0; indexJ != 3; ++indexJ)
                    {
                        UInt matIndexI = 3 * i + indexI;
                        UInt matIndexJ = 3 * j + indexJ;

                        Vec3& perturbVertexI = this->contactMesh.vertices[i];
                        Vec3& perturbVertexJ = this->contactMesh.vertices[j];

                        // Up, f(... x_i + h ... x_j + h ...)
                        perturbVertexI[indexI] += h;
                        perturbVertexJ[indexJ] += h;
                        Float upSsVal = BarrierVal();

                        // Left, f(... x_i - h ... x_j + h ...)
                        perturbVertexI[indexI] -= 2 * h;
                        Float leftSsVal = BarrierVal();

                        // Right, f(... x_i + h ... x_j - h ...)
                        perturbVertexI[indexI] += 2 * h;
                        perturbVertexJ[indexJ] -= 2 * h;
                        Float rightSsVal = BarrierVal();

                        // Down, f(... x_i - h ... x_j - h ...)
                        perturbVertexI[indexI] -= 2 * h;
                        Float downSsVal = BarrierVal();

                        Float hessEntry = (upSsVal - rightSsVal - leftSsVal + downSsVal) / (4 * hSq);
                        hess(matIndexI, matIndexJ) = hessEntry;
                    }
                }
            }
        }

        return hess;
    }

    DMat Solver::NumericalBarrierGrad_BS()
    {
        UInt size = this->bsVertCntSum;

        auto BarrierVal = [&]() -> Float
        {
            this->UpdateContactMesh();
            IPC::BuildCollisionSets(this->contactMesh, this->spatialHash, this->ground, true);
            Float val = IPC::ComputeBarrierEnergyVal(this->contactMesh, this->ground, this->contactStiffness);
            return val;
        };

        DMat grad = DMat::Zero(3 * size, 1);

        Float h = 1e-8;

        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            BSTargetInfo& info = this->targets[idx];
            for (UInt i = 0; i != size; ++i)
            {
                for (UInt indexI = 0; indexI != 3; ++indexI)
                {
                    Vec3& perturbVertex = info.target->GetControlPoint(i);
                    perturbVertex[indexI] += h;
                    Float upVal = BarrierVal();

                    perturbVertex[indexI] -= 2 * h;
                    Float downVal = BarrierVal();

                    perturbVertex[indexI] += h;
                    grad(3 * i + indexI, 0) = (upVal - downVal) / (2 * h);
                }
            }
        }

        return grad;
    }

    DMat Solver::NumericalBarrierHess_BS()
    {
        UInt size = this->bsVertCntSum;

        auto BarrierVal = [&]() -> Float
        {
            this->UpdateContactMesh();
            IPC::BuildCollisionSets(this->contactMesh, this->spatialHash, this->ground, true);
            Float val = IPC::ComputeBarrierEnergyVal(this->contactMesh, this->ground, this->contactStiffness);
            return val;
        };

        DMat hess = DMat::Zero(3 * size, 3 * size);

        Float h = 1e-8;
        Float hSq = h * h;

        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            BSTargetInfo& info = this->targets[idx];
            for (UInt i = 0; i != size; ++i)
            {
                for (UInt j = 0; j != size; ++j)
                {
                    for (UInt indexI = 0; indexI != 3; ++indexI)
                    {
                        for (UInt indexJ = 0; indexJ != 3; ++indexJ)
                        {
                            UInt matIndexI = 3 * i + indexI;
                            UInt matIndexJ = 3 * j + indexJ;

                            Vec3& perturbVertexI = info.target->GetControlPoint(i);
                            Vec3& perturbVertexJ = info.target->GetControlPoint(j);

                            // Up, f(... x_i + h ... x_j + h ...)
                            perturbVertexI[indexI] += h;
                            perturbVertexJ[indexJ] += h;
                            Float upSsVal = BarrierVal();

                            // Left, f(... x_i - h ... x_j + h ...)
                            perturbVertexI[indexI] -= 2 * h;
                            Float leftSsVal = BarrierVal();

                            // Right, f(... x_i + h ... x_j - h ...)
                            perturbVertexI[indexI] += 2 * h;
                            perturbVertexJ[indexJ] -= 2 * h;
                            Float rightSsVal = BarrierVal();

                            // Down, f(... x_i - h ... x_j - h ...)
                            perturbVertexI[indexI] -= 2 * h;
                            Float downSsVal = BarrierVal();

                            Float hessEntry = (upSsVal - rightSsVal - leftSsVal + downSsVal) / (4 * hSq);
                            hess(matIndexI, matIndexJ) = hessEntry;
                        }
                    }
                }
            }
        }

        return hess;
    }
#endif

#if defined BSIPC_NUMERIC_TEST
    //DMat Solver::NumericalEnvBarrierGrad() const
    //{
    //    UInt size = this->target->GetControlPointCnt();
    //    DMat grad = DMat::Zero(3 * size, 1);
    //    const Float h = 1e-8;

    //    for (UInt i = 0; i != 3 * size; ++i)
    //    {
    //        UInt index = i / 3;
    //        UInt coord = i % 3;

    //        Vec3& perturbVertex = this->target->GetControlPoint(index);
    //        perturbVertex[coord] += h;
    //        this->target->UpdateQuadPointCache();
    //        Float upVal = EnvBarrierVal();

    //        perturbVertex[coord] -= 2 * h;
    //        this->target->UpdateQuadPointCache();
    //        Float downVal = EnvBarrierVal();

    //        perturbVertex[coord] += h;
    //        grad(i, 0) = (upVal - downVal) / (2 * h);
    //    }

    //    return grad;
    //}

    //DMat Solver::NumericalEnvBarrierHess() const
    //{
    //    UInt size = this->target->GetControlPointCnt();
    //    const Float h = 1e-8;
    //    const Float hSq = h * h;
    //    DMat hess = DMat::Zero(3 * size, 3 * size);

    //    for (UInt i = 0; i != size; ++i)
    //    {
    //        for (UInt j = 0; j != size; ++j)
    //        {
    //            for (UInt indexI = 0; indexI != 3; ++indexI)
    //            {
    //                for (UInt indexJ = 0; indexJ != 3; ++indexJ)
    //                {
    //                    UInt matIndexI = 3 * i + indexI;
    //                    UInt matIndexJ = 3 * j + indexJ;

    //                    Vec3& perturbVertexI = this->target->GetControlPoint(i);
    //                    Vec3& perturbVertexJ = this->target->GetControlPoint(j);

    //                    // Up, f(... x_i + h ... x_j + h ...)
    //                    perturbVertexI[indexI] += h;
    //                    perturbVertexJ[indexJ] += h;
    //                    this->target->UpdateQuadPointCache();
    //                    Float upSsVal = EnvBarrierVal();

    //                    // Left, f(... x_i - h ... x_j + h ...)
    //                    perturbVertexI[indexI] -= 2 * h;
    //                    this->target->UpdateQuadPointCache();
    //                    Float leftSsVal = EnvBarrierVal();

    //                    // Right, f(... x_i + h ... x_j - h ...)
    //                    perturbVertexI[indexI] += 2 * h;
    //                    perturbVertexJ[indexJ] -= 2 * h;
    //                    this->target->UpdateQuadPointCache();
    //                    Float rightSsVal = EnvBarrierVal();

    //                    // Down, f(... x_i - h ... x_j - h ...)
    //                    perturbVertexI[indexI] -= 2 * h;
    //                    this->target->UpdateQuadPointCache();
    //                    Float downSsVal = EnvBarrierVal();

    //                    Float hessEntry = (upSsVal - rightSsVal - leftSsVal + downSsVal) / (4 * hSq);
    //                    hess(matIndexI, matIndexJ) = hessEntry;
    //                }
    //            }
    //        }
    //    }

    //    return hess;
    //}
#endif

    void Solver::FillGlobalGradient(const DMat& localGrad, BSIPC_OUT DMat& globalGrad, const std::vector<UInt> indices, tbb::mutex& mutex) const
    {
        mutex.lock();
        for (UInt i = 0; i != indices.size(); ++i)
            for (UInt coord = 0; coord != 3; ++coord)
                globalGrad(3 * indices[i] + coord, 0) += localGrad(3 * i + coord, 0);
        mutex.unlock();
    }

    void Solver::FillElasticityHessTimestepSq(
        const DMat& localHess, BSIPC_OUT SpMatData& globalHessData, const std::vector<UInt> indices, UInt patchIndex
    ) const
    {
        const Float timestepSq = this->config.timestep * this->config.timestep;

        for (UInt i = 0; i != indices.size(); ++i)
            for (UInt j = 0; j != indices.size(); ++j)
            {
                UInt globalI = indices[i], globalJ = indices[j];

                for (UInt coordI = 0; coordI != 3; ++coordI)
                    for (UInt coordJ = 0; coordJ != 3; ++coordJ)
                    {
                        // rowIndex and columnIndex are inversed cf. matrix entry indices: (i, j) <-> (columnIndex, rowIndex)
                        // row and column indices in the Hessian
                        UInt rowIndex = 3 * globalJ + coordJ;
                        UInt colIndex = 3 * globalI + coordI;

                        Float curEntryVal = timestepSq * localHess(3 * i + coordI, 3 * j + coordJ);

                        UInt inBlockIndex = (3 * i + coordI) * 27 + 3 * j + coordJ;
                        UInt globalIndex = patchIndex * 729 + inBlockIndex;

                        // Elasticity entries appear at the beginning of global Hessian, so offset = 0
                        if (this->ShouldFixHessEntry(rowIndex, colIndex))
                            globalHessData[globalIndex] = SpMatEntry(rowIndex, colIndex, 0);
                        else
                            globalHessData[globalIndex] = SpMatEntry(rowIndex, colIndex, curEntryVal);
                    }
            }
    }

    void Solver::FillSeamHess(const DMat& localHess, SpMatData& globalHessData, 
        const std::vector<UInt> localIndicesI, const std::vector<UInt> localIndicesJ, UInt bsOffsetI, UInt bsOffsetJ, 
        UInt hessEntryOffset) const
    {
        std::vector<UInt> globalIndices;
        for (UInt i = 0; i != localIndicesI.size(); ++i)
            globalIndices.emplace_back(localIndicesI[i] + bsOffsetI);
        for (UInt j = 0; j != localIndicesJ.size(); ++j)
            globalIndices.emplace_back(localIndicesJ[j] + bsOffsetJ);

        for (UInt i = 0; i != globalIndices.size(); ++i)
            for (UInt j = 0; j != globalIndices.size(); ++j)
                for (UInt coordI = 0; coordI != 3; ++coordI)
                    for (UInt coordJ = 0; coordJ != 3; ++coordJ)
                    {
                        UInt rowIndex = 3 * globalIndices[j] + coordJ;
                        UInt colIndex = 3 * globalIndices[i] + coordI;

                        Float curEntryVal = localHess(3 * i + coordI, 3 * j + coordJ);

                        UInt inBlockIndex = (3 * i + coordI) * 54 + 3 * j + coordJ;
                        UInt globalIndex = hessEntryOffset + inBlockIndex;

                        if (this->ShouldFixHessEntry(rowIndex, colIndex))
                            globalHessData[globalIndex] = SpMatEntry(rowIndex, colIndex, 0);
                        else
                            globalHessData[globalIndex] = SpMatEntry(rowIndex, colIndex, curEntryVal);
                    }
    }

    void Solver::FillLinearHess(const tbb::concurrent_vector<SpMatEntry>& trigHessData, BSIPC_OUT SpMatData& globalHessData, UInt offset) const
    {
        UInt bsHessEntryCnt = static_cast<UInt>(trigHessData.size());
        Int iterSize = 81 * bsHessEntryCnt;

        tbb::parallel_for(
            0, iterSize, 1, [&](UInt i)
            {
                UInt entryIndex = i / 81;
                UInt sourceIndexI = (i % 81) / 9;
                UInt sourceIndexJ = (i % 81) % 9;

                const SpMatEntry& entry = trigHessData[entryIndex];
                UInt vertexIndexI = entry.row() / 3;
                UInt vertexIndexJ = entry.col() / 3;
                UInt vertexCompI = entry.row() % 3;
                UInt vertexCompJ = entry.col() % 3;

                const SourceInfo& sourceI = this->linearVertexSources[vertexIndexI][sourceIndexI];
                const SourceInfo& sourceJ = this->linearVertexSources[vertexIndexJ][sourceIndexJ];

                Float coeffI = sourceI.coeff, coeffJ = sourceJ.coeff;

                UInt targetIndexI = sourceI.index, targetIndexJ = sourceJ.index;

                UInt rowIndex = 3 * targetIndexI + vertexCompI;
                UInt colIndex = 3 * targetIndexJ + vertexCompJ;
                Float val = entry.value() * coeffI * coeffJ;

                if (this->ShouldFixHessEntry(rowIndex, colIndex))
                    globalHessData[offset + i] = SpMatEntry(rowIndex, colIndex, 0);
                else
                    globalHessData[offset + i] = SpMatEntry(rowIndex, colIndex, val);
            }
        );
    }

    void Solver::FillAnimatedLinearHess(
        const tbb::concurrent_vector<SpMatEntry>& bsAnimHessData, const tbb::concurrent_vector<SpMatEntry>& animBsHessData, 
        const tbb::concurrent_vector<SpMatEntry>& animAnimHessData, SpMatData& globalHessData, UInt offset) const
    {
        UInt bsAnimHessEntryCnt = 9 * bsAnimHessData.size();                // rows are bs vertices
        UInt animBsHessEntryCnt = 9 * animBsHessData.size();                // cols are bs vertices
        UInt animAnimHessEntryCnt = animAnimHessData.size();

        tbb::parallel_for(
            0, static_cast<Int>(bsAnimHessEntryCnt), 1, [&](UInt i)
        //for (UInt i = 0; i != bsAnimHessEntryCnt; ++i)
            {
                UInt entryIndex = i / 9;
                UInt sourceIndex = i % 9;

                const SpMatEntry& entry = bsAnimHessData[entryIndex];
                UInt vertexIndex = entry.row() / 3;
                UInt vertexComp = entry.row() % 3;

                const SourceInfo& source = this->linearVertexSources[vertexIndex][sourceIndex];
                Float coeff = source.coeff;
                UInt targetIndex = source.index;

                UInt rowIndex = 3 * targetIndex + vertexComp;
                UInt colIndex = entry.col() - 3 * this->trigVertCntSum + 3 * this->bsVertCntSum;
                Float val = entry.value() * coeff;

                if (this->ShouldFixHessEntry(rowIndex, colIndex))
                    globalHessData[offset + i] = SpMatEntry(rowIndex, colIndex, 0);
                else
                    globalHessData[offset + i] = SpMatEntry(rowIndex, colIndex, val);
            }
        );

        offset += bsAnimHessEntryCnt;

        tbb::parallel_for(
            0, static_cast<Int>(animBsHessEntryCnt), 1, [&](UInt i)
        //for (UInt i = 0; i != animBsHessEntryCnt; ++i)
            {
                UInt entryIndex = i / 9;
                UInt sourceIndex = i % 9;

                const SpMatEntry& entry = animBsHessData[entryIndex];
                UInt vertexIndex = entry.col() / 3;
                UInt vertexComp = entry.col() % 3;

                const SourceInfo& source = this->linearVertexSources[vertexIndex][sourceIndex];
                Float coeff = source.coeff;
                UInt targetIndex = source.index;

                UInt rowIndex = entry.row() - 3 * this->trigVertCntSum + 3 * this->bsVertCntSum;
                UInt colIndex = 3 * targetIndex + vertexComp;
                Float val = entry.value() * coeff;

                if (this->ShouldFixHessEntry(rowIndex, colIndex))
                    globalHessData[offset + i] = SpMatEntry(rowIndex, colIndex, 0);
                else
                    globalHessData[offset + i] = SpMatEntry(rowIndex, colIndex, val);
            }
        );

        offset += animBsHessEntryCnt;

        tbb::parallel_for(
            0, static_cast<Int>(animAnimHessEntryCnt), 1, [&](UInt i)
            {
                const SpMatEntry& entry = animAnimHessData[i];
                UInt rowIndex = entry.row() - 3 * this->trigVertCntSum + 3 * this->bsVertCntSum;
                UInt colIndex = entry.col() - 3 * this->trigVertCntSum + 3 * this->bsVertCntSum;

                globalHessData[offset + i] = SpMatEntry(rowIndex, colIndex, entry.value());
            }
        );
    }

    inline void Solver::DBCOnGradient(DMat& gradient)
    {
        BSIPC_ASSERT(
            gradient.cols() == 1 && gradient.rows() == this->bsVertCntSum * 3, 
            "Gradient size do not match"
        );

        for (const BSTargetInfo& info : this->targets)
        {
            UInt offset = info.bsVertIndexOffset.value();
            const std::unordered_set<UInt>& indices = this->config.surfacesInfo[info.index].fixedNodes;
            for (UInt index : indices)
                for (UInt coord = 0; coord != 3; ++coord)
                    gradient(3 * (index + offset) + coord, 0) = 0;

            const std::vector<MovingNode>& movingNodes = this->config.surfacesInfo[info.index].movingNodes;
            for (const MovingNode& node : movingNodes)
            {
                UInt index = node.index;
                for (UInt coord = 0; coord != 3; ++coord)
                    gradient(3 * (index + offset) + coord, 0) = 0;
            }
        }
    }

    Bool Solver::ShouldFixHessEntry(UInt indexI, UInt indexJ) const
    {
        UInt size = this->bsVertCntSum;
        const std::vector<BSTargetInfo>& infos = this->targets;

        //BSIPC_ASSERT(
        //    indexI < size * 3 && indexJ < size * 3,
        //    "Indices out of bounds"
        //);

        std::vector<UInt> bsStartIndex(infos.size());
        for (UInt i = 0; i < infos.size(); i++)
            bsStartIndex[i] = infos[i].bsVertIndexOffset.value();

        UInt surfIndexI = 0, surfIndexJ = 0;

        for (UInt i = 0; i != bsStartIndex.size(); ++i)
        {
            if (indexI >= bsStartIndex[i] * 3)
                surfIndexI = i;
            if (indexJ >= bsStartIndex[i] * 3)
                surfIndexJ = i;
        }

        const BSSurfaceInfo& surfInfoI = config.surfacesInfo[surfIndexI];
        const BSSurfaceInfo& surfInfoJ = config.surfacesInfo[surfIndexJ];

        UInt bsVertOffsetI = infos[surfIndexI].bsVertIndexOffset.value();
        UInt bsVertOffsetJ = infos[surfIndexJ].bsVertIndexOffset.value();
        UInt localIndexI = indexI - bsVertOffsetI * 3;
        UInt localIndexJ = indexJ - bsVertOffsetJ * 3;

        UInt ctrlIndexI = localIndexI / 3;
        UInt ctrlIndexJ = localIndexJ / 3;

        return surfInfoI.fixedNodes.contains(ctrlIndexI) || surfInfoJ.fixedNodes.contains(ctrlIndexJ)
            || (config.MovingNodeVelocity(surfIndexI, ctrlIndexI) != std::nullopt)
            || (config.MovingNodeVelocity(surfIndexJ, ctrlIndexJ) != std::nullopt);
    }

    Bool Solver::ActivateAnimatedMesh() const
    {
        Bool result = false;
        if (this->animatedMeshCache.has_value() && !this->animatedMeshPoses.empty())
        {
            if (this->stepCache.CurStep() < this->config.animatedMeshInfo.value().length)
            {
                UInt curVertCnt = this->animatedMeshCache.value().vertCnt;
                DVec curPos = DVec::Zero(3 * curVertCnt);
                DVec targetPos = DVec::Zero(3 * curVertCnt);

                tbb::parallel_for(
                    0, static_cast<Int>(curVertCnt), 1, [&](UInt i)
                    {
                        Vec3 curVert = this->contactMesh.vertices[this->animatedMeshCache.value().startIdx + i];
                        Vec3 targetVert = this->animatedMeshPoses[this->stepCache.CurStep()].vertices[i];
                        for (UInt coord = 0; coord != 3; ++coord)
                        {
                            curPos(3 * i + coord) = curVert[coord];
                            targetPos(3 * i + coord) = targetVert[coord];
                        }
                    }
                );

                DVec diff = curPos - targetPos;
                if (diff.lpNorm<Eigen::Infinity>() > this->config.animatedMeshInfo.value().tolerance)
                    result = true;

                this->AppendToLog(fmt::format("diff.inf_norm: {:.6e}, diff.l2_norm: {:.6e}. ", diff.lpNorm<Eigen::Infinity>(), diff.norm()));
            }
        }
        return result;
    }

    Float Solver::IPVal(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos)
    {
        Float accPotEnergyVal = 0;
        for (const BSTargetInfo& info : this->targets)
        {
            UInt ssIterSize = static_cast<UInt>(info.target->GetSsQuadPoints().size());
            Float ssEnergyVal =
                oneapi::tbb::parallel_reduce
                (
                    oneapi::tbb::blocked_range<UInt>(0, ssIterSize),
                    0.,
                    [&](const oneapi::tbb::blocked_range<UInt>& range, Float localSum) -> Float
                    {
                        for (UInt i = range.begin(); i != range.end(); ++i)
                        {
                            const QuadPoint& curQuadPoint = info.target->GetSsQuadPoints()[i];

                            Float weight = curQuadPoint.Weight();

                            Float ssDensity = Energy::SsValAt(curQuadPoint, this->config, info.target);
                            //Float ssDensity = Energy::StVKSsValAt(curQuadPoint, this->config, info.target);
                            //Float bdDensity = Energy::BdValAt(curQuadPoint, this->config, info.target);

                            localSum += weight * (ssDensity);
                            //localSum += weight * (bdDensity);
                        }

                        return localSum;
                    },
                    [&](Float a, Float b) -> Float
                    {
                        return a + b;
                    }
                );
            this->AppendToLog(fmt::format("ssEnergyVal: {:.6e}. ", ssEnergyVal));

            UInt bdIterSize = static_cast<UInt>(info.target->GetBdQuadPoints().size());
            Float bdEnergyVal =
                oneapi::tbb::parallel_reduce
                (
                    oneapi::tbb::blocked_range<UInt>(0, bdIterSize),
                    0.,
                    [&](const oneapi::tbb::blocked_range<UInt>& range, Float localSum) -> Float
                    {
                        for (UInt i = range.begin(); i != range.end(); ++i)
                        {
                            const QuadPoint& curQuadPoint = info.target->GetBdQuadPoints()[i];

                            Float weight = curQuadPoint.Weight();

                            Float bdDensity = Energy::BdValAt(curQuadPoint, this->config, info.target);

                            localSum += weight * (bdDensity);
                        }

                        return localSum;
                    },
                    [&](Float a, Float b) -> Float
                    {
                        return a + b;
                    }
                );

            accPotEnergyVal += ssEnergyVal;
            accPotEnergyVal += bdEnergyVal;
        }

        Float inertiaVal = Energy::InertiaVal(hypPos, estPos, this->massMatDiagEntries, this->config.thickness);
        Float barrierVal = IPC::ComputeBarrierEnergyVal(this->contactMesh, this->ground, this->contactStiffness);
        Float frictionVal = IPC::ComputeFrictionEnergyVal(this->contactMesh, this->ground);

        Float timestepSquared = this->config.timestep * this->config.timestep;

        Float IPVal = inertiaVal + barrierVal + frictionVal + timestepSquared * accPotEnergyVal;
        //Float IPVal = barrierVal + frictionVal + timestepSquared * accPotEnergyVal;

        if (this->includeAnimatedMeshInSystem)
        {
            Float kappa = this->config.animatedMeshInfo.value().stiffness;
            Float penalty = 0;
            const std::vector<Vec3>& nextFramePos = this->animatedMeshPoses[this->stepCache.CurStep()].vertices;
            const std::vector<Vec3>& curPos = this->animatedMeshCache.value().curPos;

            BSIPC_ASSERT(nextFramePos.size() == curPos.size(), "animated mesh info has inconsistent vertex size");
            for (UInt i = 0; i != nextFramePos.size(); ++i)
            {
                Vec3 diff = curPos[i] - nextFramePos[i];
                penalty += 0.5 * kappa * diff.squaredNorm();
            }

            IPVal += penalty;
        }

        return IPVal;
    }

    DMat Solver::IPGrad(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos)
    {
        UInt size = this->bsVertCntSum;
        DMat localElasticityGrad = DMat::Zero(3 * size, 1);
        UInt surfPatchStartIndex = 0;

        for (const BSTargetInfo& info : this->targets)
        {
            Int ssIterSize = static_cast<Int>(info.target->GetSsQuadPoints().size());

            tbb::mutex mutex;

            UInt patchSize = info.target->GetPatchCnt();
            const std::vector<UInt>& patchStartIndices = info.target->GetPatchStartIndices();
            const std::vector<QuadPoint>& curSsQuadPoints = info.target->GetSsQuadPoints();
            UInt curOffset = info.bsVertIndexOffset.value();

            tbb::parallel_for(
                0, static_cast<Int>(patchSize), 1, [&](Int idx)
                //for (Int idx = 0; idx != patchSize; ++idx)
                {
                    UInt patchQuadStartIndex = patchStartIndices[idx];

                    UInt patchQuadCnt = 0;
                    if (idx == patchStartIndices.size() - 1)
                        patchQuadCnt = info.target->GetSsQuadPoints().size() - patchQuadStartIndex;
                    else
                        patchQuadCnt = patchStartIndices[idx + 1] - patchQuadStartIndex;

                    Mat<27, 1> patchBlockGrad = Mat<27, 1>::Zero();

                    const QuadPoint& startQuadPoint = curSsQuadPoints[patchQuadStartIndex];
                    std::vector<UInt> localIndices = info.target->SupportedNodesAt(info.target->GetBdQuadPoints()[idx]);

                    std::vector<UInt> globalIndices(localIndices.size());
                    for (UInt i = 0; i != localIndices.size(); ++i)
                        globalIndices[i] = localIndices[i] + curOffset;

                    // membrane gradient
                    for (UInt i = 0; i != patchQuadCnt; ++i)
                    {
                        const QuadPoint& curSsQuadPoint = curSsQuadPoints[patchQuadStartIndex + i];

                        if (curSsQuadPoint.U() == 0 && curSsQuadPoint.V() == 0)
                            continue;

                        Float weight = curSsQuadPoint.Weight();

                        Mat<27, 1> ssBlockGrad = weight * Energy::SsGradBlockAt(curSsQuadPoint, localIndices, this->config, info.target);
                        //Mat<27, 1> ssBlockGrad = weight * Energy::StVKSsGradBlockAt(curSsQuadPoint, localIndices, this->config, info.target);
                        patchBlockGrad += ssBlockGrad;
                    }

                    // bending gradient
                    const QuadPoint& curBdQuadPoint = info.target->GetBdQuadPoints()[idx];
                    {
                        Float weight = curBdQuadPoint.Weight();
                        Mat<27, 1> bdBlockGrad = weight * Energy::BdGradBlockAt(curBdQuadPoint, localIndices, this->config, info.target);
                        patchBlockGrad += bdBlockGrad;
                    }

                    this->localElasticityGrad[surfPatchStartIndex + idx] = patchBlockGrad;
                    this->FillGlobalGradient(patchBlockGrad, localElasticityGrad, globalIndices, mutex);
                }
            );

            surfPatchStartIndex += patchSize;
        }

#if defined BSIPC_NUMERIC_TEST
        //DMat numericalSeamGrad = Energy::NumericalSeamGrad(this->config, this->targets);

        //DMat seamGrad_T = seamGrad.transpose();
        //DMat numericalSeamGrad_T = numericalSeamGrad.transpose();
        //DMat diff = seamGrad_T - numericalSeamGrad_T;

        //this->AppendToLog("===== Seam Grad Test (Global) =====\n\n");
        //this->AppendToLog("Analytic Seam Grad: \n" + ToStr(seamGrad_T) + "\n");
        //this->AppendToLog("Numerical Seam Grad: \n" + ToStr(numericalSeamGrad_T) + "\n");
        //this->AppendToLog("Diff: \n" + ToStr(diff) + "\n");
        //this->AppendToLog(fmt::format("diff.norm_inf = {}", diff.lpNorm<Eigen::Infinity>()));
        //this->AppendToLog("\n<<<<< End Test >>>>>\n");
#endif
        DMat inertiaGrad = Energy::InertiaGrad(hypPos, estPos, this->massMatDiagEntries, this->config.thickness);
        DMat barrierGrad = this->AccumulateTrigGrad(this->contactMeshBarrierGrad);

        Float timestepSquared = this->config.timestep * this->config.timestep;
        DMat IPGrad = inertiaGrad + barrierGrad + timestepSquared * localElasticityGrad;
        //DMat IPGrad = barrierGrad + timestepSquared * localElasticityGrad;

        DBCOnGradient(IPGrad);

#if defined BSIPC_NUMERIC_TEST
        //DMat numericalIPGrad = this->NumericalIPGrad(hypPos, estPos);

        //DMat IPGrad_T = IPGrad.transpose();
        //DMat numericalIPGrad_T = numericalIPGrad.transpose();

        //this->AppendToLog("===== IPGrad Test =====\n\n");
        //this->AppendToLog("Analytic IPGrad: \n" + ToStr(IPGrad_T) + "\n");
        //this->AppendToLog("Numerical IPGrad: \n" + ToStr(numericalIPGrad_T) + "\n");
        //this->AppendToLog("\n<<<<< End Test >>>>>\n");

        //DMat numericalSeamGrad = Energy::NumericalSeamGrad(this->config, this->targets);

        //DMat seamGrad_T = seamGrad.transpose();
        //DMat numericalSeamGrad_T = numericalSeamGrad.transpose();
        //DMat diff = seamGrad_T - numericalSeamGrad_T;

        //for (UInt i = 0; i != diff.rows(); ++i)
        //    if (std::abs(diff(i, 0)) < 1e-5)
        //        diff(i, 0) = 0;

        //this->AppendToLog("===== Seam Grad Test (Assembled) =====\n\n");
        //this->AppendToLog("Analytic Seam Grad: \n" + ToStr(seamGrad_T) + "\n");
        //this->AppendToLog("Numerical Seam Grad: \n" + ToStr(numericalSeamGrad_T) + "\n");
        //this->AppendToLog("Diff: \n" + ToStr(diff) + "\n");
        //this->AppendToLog("\n<<<<< End Test >>>>>\n");
#endif

        if (this->includeAnimatedMeshInSystem)
        {
            DMat amalgGrad = DMat::Zero(IPGrad.size() + 3 * this->animatedMeshCache.value().vertCnt, 1);

            // preserve previous IPGrad
            UInt offset = IPGrad.rows();
            amalgGrad.topRows(IPGrad.rows()) = IPGrad;

            // include penalty
            const std::vector<Vec3>& nextFramePos = this->animatedMeshPoses[this->stepCache.CurStep()].vertices;
            const std::vector<Vec3>& curPos = this->animatedMeshCache.value().curPos;
            for (UInt i = 0; i != nextFramePos.size(); ++i)
            {
                Vec3 diff = curPos[i] - nextFramePos[i];
                for (UInt coord = 0; coord != 3; ++coord)
                    amalgGrad(offset + 3 * i + coord, 0) = this->config.animatedMeshInfo.value().stiffness * diff[coord];
            }

            // include IPC terms
            tbb::parallel_for(
                0, static_cast<Int>(this->animatedMeshCache.value().vertCnt), 1, [&](Int i)
                {
                    amalgGrad.block<3, 1>(offset + 3 * i, 0) += this->contactMeshBarrierGrad[this->animatedMeshCache.value().startIdx + i];
                }
            );

            IPGrad = std::move(amalgGrad);
        }

        return IPGrad;
    }

    SpMat Solver::IPHess(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos)
    {
        //SpMat _contactBarrierHess;
        //BSIPC_PROFILE_START("Optimized Barrier Hess Conversion (w/o spatial encoding)");
        //this->_ContactBarrierHess(_contactBarrierHess);
        //BSIPC_PROFILE_END("Optimized Barrier Hess Conversion (w/o spatial encoding)");

        SpMat contactBarrierHess;
        BSIPC_PROFILE_START("Optimized Barrier Hess Conversion");
        this->ContactBarrierHessReordered(contactBarrierHess);
        //ValidateSparseMatrix(contactBarrierHess, true);
        BSIPC_PROFILE_END("Optimized Barrier Hess Conversion");

        //SpMat _hess;
        //BSIPC_PROFILE_START("Optimized Barrier Hess Conversion Opt");
        //this->ContactBarrierHessReorderedOpt(_hess);
        //BSIPC_PROFILE_END("Optimized Barrier Hess Conversion Opt");

        //SpMat diff = contactBarrierHess - _hess;
        //Float maxDiff = 0.;
        //for (Int k = 0; k != diff.outerSize(); ++k)
            //for (SpMat::InnerIterator it(diff, k); it; ++it)
            //{
                //if (it.value() > 1e-5)
                    //BSIPC_PEEK(it.value());
                //maxDiff = std::max(maxDiff, std::abs(it.value()));
            //}
        //BSIPC_PEEK(maxDiff);

        //SpMat _1;
        //BSIPC_PROFILE_START("Barrier Hess Conversion Triplet");
        //this->_ContactBarrierHessTriplet(contactBarrierHess);
        //BSIPC_PROFILE_END("Barrier Hess Conversion Triplet");

        BSIPC_PROFILE_START("Hess Prep");

        // No need to take into consideration the gravity term as its Hessian is zero
        UInt size = this->bsVertCntSum;

        UInt trigMeshVertCnt = 0, quadPointsCnt = 0, patchCnt = 0;
        for (const BSTargetInfo& info : this->targets)
        {
            quadPointsCnt += info.target->GetSsQuadPoints().size();
            trigMeshVertCnt += info.trigVertCnt.value();
            patchCnt += info.target->UDomainMax() * info.target->VDomainMax();
        }

        UInt seamPointPairCnt = 0;
        for (const SeamInfo& info : this->config.seamInfos)
            seamPointPairCnt += info.seamPoints.size();

        BSIPC_PROFILE_START("Hess Data Allocation");

        // Layout:
        //  Potential Energy Hessian | Inertia Hessian | Seam Penalty Hessian

        UInt hessEntryCnt = this->inertiaHessEntryCnt + this->potHessEntryCnt;
        if (this->includeAnimatedMeshInSystem)
            hessEntryCnt += 3 * this->animatedMeshCache.value().vertCnt;

        this->hessEntries.resize(hessEntryCnt);

        tbb::parallel_for(
            0, static_cast<Int>(this->hessEntries.size()), 1, [&](Int i)
            {
                this->hessEntries[i] = SpMatEntry(0, 0, 0);
            }
        );

        BSIPC_PROFILE_END("Hess Data Allocation");

        BSIPC_PROFILE_END("Hess Prep");

        BSIPC_PROFILE_START("Elasticity Hessian Blocks");

        UInt quadPtStartIndex = 0;
        UInt surfPatchStartIndex = 0;

        for (const BSTargetInfo& info : this->targets)
        {
            const BSSurface& target = *info.target;
            const std::vector<QuadPoint>& curSsQuadPoints = target.GetSsQuadPoints();

            UInt curOffset = info.bsVertIndexOffset.value();
            const Float timestepSq = this->config.timestep * this->config.timestep;
            UInt patchSize = target.GetPatchCnt();

#if defined BSIPC_NUMERIC_TEST
            //for (UInt i = 0; i != curQuadPoints.size(); ++i)
            //{
            //    const QuadPoint& curQuadPoint = info.target->GetQuadPoints()[i];

            //    Float weight = curQuadPoint.Weight();

            //    std::vector<UInt> localIndices = info.target->SupportedNotesAt(curQuadPoint);
            //    DMat ssBlockHess = weight * Energy::SsHessBlockAt(curQuadPoint, localIndices, this->config, info.target);
            //    DMat bdBlockHess = weight * Energy::BdHessBlockAt(curQuadPoint, localIndices, this->config, info.target);

            //    this->SsHessBlockNumericTest(curQuadPoint, localIndices, Energy::SsHessBlockAt(curQuadPoint, localIndices, this->config, info.target), info);

            //    DMat localBlockHess = ssBlockHess + bdBlockHess;
            //    //DMat localBlockHess = bdBlockHess;

            //    UInt curOffset = info.bsVertIndexOffset.value();

            //    std::vector<UInt> globalIndices(localIndices.size());
            //    for (UInt idx = 0; idx != localIndices.size(); ++idx)
            //        globalIndices[idx] = localIndices[idx] + curOffset;

            //    FillElasticityHessTimestepSq(localBlockHess, hessData, 0, globalIndices, i, info);
            //    //FillElasticityHessTimestepSq(bdBlockHess, bdData, 0, indices, i);
            //}
#endif
            const std::vector<UInt>& patchStartIndices = target.GetPatchStartIndices();
            BSIPC_ASSERT(
                patchSize == patchStartIndices.size(), "Size of [patchStartIndices] should match the number of patches"
            );

            BSIPC_PEEK(curSsQuadPoints.size());

            tbb::parallel_for(
                0, static_cast<Int>(patchSize), 1, [&](Int idx)
                //for (Int idx = 0; idx != patchSize; ++idx)
                {
                    UInt patchQuadStartIndex = patchStartIndices[idx];

                    UInt patchQuadCnt = 0;
                    if (idx == patchStartIndices.size() - 1)
                        patchQuadCnt = target.GetSsQuadPoints().size() - patchQuadStartIndex;
                    else
                        patchQuadCnt = patchStartIndices[idx + 1] - patchQuadStartIndex;

                    Mat<27, 27> patchBlockHess = Mat<27, 27>::Zero();

                    const QuadPoint& startQuadPoint = curSsQuadPoints[patchQuadStartIndex];
                    std::vector<UInt> localIndices = info.target->SupportedNodesAt(info.target->GetBdQuadPoints()[idx]);

                    std::vector<UInt> globalIndices(localIndices.size());
                    for (UInt i = 0; i != localIndices.size(); ++i)
                        globalIndices[i] = localIndices[i] + curOffset;

                    for (UInt i = 0; i != patchQuadCnt; ++i)
                    {
                        const QuadPoint& curSsQuadPoint = curSsQuadPoints[patchQuadStartIndex + i];

                        if (curSsQuadPoint.U() == 0 && curSsQuadPoint.V() == 0)
                            continue;

                        Float weight = curSsQuadPoint.Weight();

                        Mat<27, 27> ssBlockHess = weight * Energy::SsHessBlockAt(curSsQuadPoint, localIndices, this->config, info.target);
                        //Mat<27, 27> ssBlockHess = weight * Energy::StVKSsHessBlockAt(curSsQuadPoint, localIndices, this->config, info.target);
                        patchBlockHess += ssBlockHess;
                    }
                    patchBlockHess += this->localBendingHess[surfPatchStartIndex + idx];
                    this->localElasticityHess[surfPatchStartIndex + idx] = patchBlockHess;

                    FillElasticityHessTimestepSq(patchBlockHess, this->hessEntries, globalIndices, surfPatchStartIndex + idx);
                }
            );

            quadPtStartIndex += curSsQuadPoints.size();
            surfPatchStartIndex += patchSize;
        }

        BSIPC_PROFILE_END("Elasticity Hessian Blocks");

        BSIPC_PROFILE_SCOPE(
            "Inertia Hess",
            Energy::FillInertiaHess(this->hessEntries, potHessEntryCnt, this->massMatDiagEntries, this->config, *this);
        )

#if defined BSIPC_NUMERIC_TEST
        //DMat numericalSeamHess = Energy::NumericalSeamHess(this->config, this->targets);
        //SpMat seamHessMat(3 * size, 3 * size);
        //seamHessMat.setFromTriplets(seamHessData.begin(), seamHessData.end());
        //DMat seamHessD = SparseToDense(seamHessMat);
        //DMat diff = seamHessD - numericalSeamHess;

        //this->AppendToLog("===== Seam Hessian Test (Assembled) =====\n\n");
        //this->AppendToLog("Analytic Seam Hessian: \n" + ToStr(seamHessD) + "\n");
        //this->AppendToLog("Numerical Seam Hessian: \n" + ToStr(numericalSeamHess) + "\n");
        //this->AppendToLog("Diff: \n" + ToStr(diff) + "\n");
        //this->AppendToLog(fmt::format("diff.norm_inf = {}", diff.lpNorm<Eigen::Infinity>()));
        //this->AppendToLog("\n<<<<< End Test >>>>>\n");
#endif

        BSIPC_PEEK(this->hessEntries.size());
        Float dt = this->config.timestep;

        /////// TEST /////

        //UInt elasInertiaHessSize = 3 * size;
        //SpMat testHess(elasInertiaHessSize, elasInertiaHessSize);
        ////testHess.setFromTriplets(this->hessEntries.begin(), this->hessEntries.begin() + this->inertiaHessEntryCnt + this->potHessEntryCnt);
        //testHess.setFromTriplets(this->hessEntries.begin(), this->hessEntries.begin() + this->inertiaHessEntryCnt);
        //testHess.makeCompressed();

        //this->elasticityHess.resize(this->elasticityHessTemplate.rows(), this->elasticityHessTemplate.cols());
        //this->elasticityHess.reserve(this->elasticityHessTemplate.nonZeros());
        //std::memcpy(this->elasticityHess.valuePtr(), this->elasticityHessTemplate.valuePtr(), this->elasticityHessTemplate.nonZeros() * sizeof(Float));
        //std::memcpy(this->elasticityHess.innerIndexPtr(), this->elasticityHessTemplate.innerIndexPtr(), this->elasticityHessTemplate.nonZeros() * sizeof(UInt));
        //std::memcpy(this->elasticityHess.outerIndexPtr(), this->elasticityHessTemplate.outerIndexPtr(), (this->elasticityHessTemplate.outerSize() + 1) * sizeof(UInt));
        //this->elasticityHess.resizeNonZeros(this->elasticityHessTemplate.nonZeros());
        //this->elasticityHess.makeCompressed();

        ////tbb::parallel_for(
        ////    0, static_cast<Int>(this->elasticityHess.nonZeros()), 1, [&](Int i)
        ////    {
        ////        for (const ElasticityHessEntrySource& source : this->elasticityHessSources[i])
        ////        {
        ////            this->elasticityHess.valuePtr()[i] += dt * dt * this->localElasticityHess[source.globalPatchIdx](source.rowIdx, source.colIdx);
        ////        }
        ////    }
        ////);

        //// get the maximum entry in testHess
        //Float maxTestEntry = 0.;
        //for (Int k = 0; k != testHess.outerSize(); ++k)
        //    for (SpMat::InnerIterator it(testHess, k); it; ++it)
        //        maxTestEntry = std::max(maxTestEntry, std::abs(it.value()));

        //Float maxElasticityEntry = 0.;
        //for (Int k = 0; k != this->elasticityHess.outerSize(); ++k)
        //    for (SpMat::InnerIterator it(this->elasticityHess, k); it; ++it)
        //        maxElasticityEntry = std::max(maxElasticityEntry, std::abs(it.value()));

        //BSIPC_INFO("testHess: {}x{}, this->elasticityHess: {}x{}", testHess.rows(), testHess.cols(), this->elasticityHess.rows(), this->elasticityHess.cols());
        //SpMat diff = testHess - this->elasticityHess;

        //this->AppendToLog("gt: \n");
        //this->AppendToLog(ToStr(SparseToDense(testHess)));
        //this->AppendToLog("new assemble: \n");
        //this->AppendToLog(ToStr(this->elasticityHessTemplate));
        ////this->AppendToLog(ToStr(SparseToDense(this->elasticityHess)));

        //Float maxEntry = 0.;
        //for (Int k = 0; k != diff.outerSize(); ++k)
        //    for (SpMat::InnerIterator it(diff, k); it; ++it)
        //        maxEntry = std::max(maxEntry, std::abs(it.value()));
        //BSIPC_WARN("max entry diff: {} // max gt entry: {} // max test entry: {}", maxEntry, maxTestEntry, maxElasticityEntry);

        /////// END TEST /////

        if (this->includeAnimatedMeshInSystem)
        {
            size += this->animatedMeshCache.value().vertCnt;

            //UInt hessOffset = this->inertiaHessEntryCnt + this->potHessEntryCnt + barrierHessEntryCnt + seamHessEntryCnt;
            UInt hessOffset = this->inertiaHessEntryCnt + this->potHessEntryCnt;
            
            tbb::parallel_for(
                0, static_cast<Int>(3 * this->animatedMeshCache.value().vertCnt), 1, [&](Int i)
                {
                    this->hessEntries[hessOffset + i] = SpMatEntry(3 * this->bsVertCntSum + i, 3 * this->bsVertCntSum + i,
                        this->config.animatedMeshInfo.value().stiffness);
                }
            );
        }

        BSIPC_PROFILE_SCOPE(
            "Hess Assembly",
            UInt hessSize = 3 * size;
            SpMat IPHess(hessSize, hessSize);
            IPHess.setFromTriplets(this->hessEntries.begin(), this->hessEntries.end());
        );

        // Test Hessian composition
#if defined BSIPC_NUMERIC_TEST
        // {TEST} when testing, remove the SPD projection in bending hessian computing
        //this->BdHessNumericTest(bdData);
#endif

#if defined BSIPC_NUMERIC_TEST
        //DMat numericalIPHess = this->NumericalIPHess(hypPos, estPos);

        //DMat IPHess_D = SparseToDense(IPHess);

        //DMat diff = IPHess_D - numericalIPHess;

        //for (UInt i = 0; i != diff.cols(); ++i)
        //    for (UInt j = 0; j != diff.cols(); ++j)
        //        if (std::abs(diff(i, j)) < 1e-5)
        //            diff(i, j) = 0;

        //this->AppendToLog("===== IPHess Test =====\n\n");
        //this->AppendToLog("Analytic IPGrad: \n" + ToStr(IPHess_D) + "\n");
        //this->AppendToLog("Numerical IPGrad: \n" + ToStr(numericalIPHess) + "\n");
        //this->AppendToLog("diff: \n" + ToStr(diff) + "\n");
        //this->AppendToLog("\n<<<<< End Test >>>>>\n");
#endif

        BSIPC_PROFILE_START("Sum all parts of Hess");
        IPHess = IPHess + contactBarrierHess;
        BSIPC_PROFILE_END("Sum all parts of Hess");

        return IPHess;
    }

    void Solver::_ContactBarrierHess(SpMat& hess)
    {
        BSIPC_PROFILE_START("[Func] ContactBarrierHessReordered");

        BSIPC_PROFILE_START("Init");
        const IPC::BHessian& barrierHess = this->contactMeshBarrierHess;
        auto CompressKey = [](UInt i, UInt j) -> uint64_t {
            return (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
        };
        auto DecompressHash = [](uint64_t key) -> std::pair<UInt, UInt> {
            UInt i = static_cast<UInt>(key >> 32);
            UInt j = static_cast<UInt>(key & 0xFFFFFFFF);
            return { i, j };
        };

        UInt size = this->bsVertCntSum;
        if (this->includeAnimatedMeshInSystem)
            size += this->animatedMeshCache.value().vertCnt;

        using LinearHessBlockCache = std::pair<Mat<3, 3>, tbb::concurrent_vector<LinearHessBlockSource>>;
        tbb::concurrent_unordered_map<uint64_t, LinearHessBlockCache> linearBlockSourceMap;

        BSIPC_PROFILE_END("Init");

        BSIPC_PROFILE_START("Record Linear Hess Block Sources");
        {
            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.D2Index.size()), 1, [&](UInt idx)
                {
                    for (UInt i = 0; i != 2; ++i)
                        for (UInt j = 0; j != 2; ++j)
                        {
                            uint64_t key = CompressKey(barrierHess.D2Index[idx][i], barrierHess.D2Index[idx][j]);
                            LinearHessBlockSource source{
                                .blockType = LinearHessBlockType::Block6x6,
                                .vecIndex = idx,
                                .blockIndices = {i, j}
                            };

                            auto insertResult = linearBlockSourceMap.insert({ 
                                key, 
                                std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{}) 
                            });

                            insertResult.first->second.second.emplace_back(source);
                        }
                }
            );

            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.D3Index.size()), 1, [&](UInt idx)
                {
                    for (UInt i = 0; i != 3; ++i)
                        for (UInt j = 0; j != 3; ++j)
                        {
                            uint64_t key = CompressKey(barrierHess.D3Index[idx][i], barrierHess.D3Index[idx][j]);
                            LinearHessBlockSource source{
                                .blockType = LinearHessBlockType::Block9x9,
                                .vecIndex = idx,
                                .blockIndices = {i, j}
                            };

                            auto insertResult = linearBlockSourceMap.insert({
                                key,
                                std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{})
                            });

                            insertResult.first->second.second.emplace_back(source);
                        }
                }
            );

            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.D4Index.size()), 1, [&](UInt idx)
                {
                    for (UInt i = 0; i != 4; ++i)
                        for (UInt j = 0; j != 4; ++j)
                        {
                            uint64_t key = CompressKey(barrierHess.D4Index[idx][i], barrierHess.D4Index[idx][j]);
                            LinearHessBlockSource source{
                                .blockType = LinearHessBlockType::Block12x12,
                                .vecIndex = idx,
                                .blockIndices = {i, j}
                            };

                            auto insertResult = linearBlockSourceMap.insert({
                                key,
                                std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{})
                            });

                            insertResult.first->second.second.emplace_back(source);
                        }
                }
            );
        }
        BSIPC_PROFILE_END("Record Linear Hess Block Sources");

        BSIPC_PROFILE_START("Accumulate Linear Hess Blocks");
        {
            tbb::parallel_for_each(
                linearBlockSourceMap.begin(), linearBlockSourceMap.end(),
                [&](auto& block)
                {
                    Mat<3, 3> accBlock = Mat<3, 3>::Zero();
                    for (const LinearHessBlockSource& source : block.second.second)
                    {
                        switch (source.blockType)
                        {
                        case LinearHessBlockType::Block3x3:
                            accBlock += barrierHess.H3x3[source.vecIndex];
                            break;
                        case LinearHessBlockType::Block6x6:
                            accBlock += barrierHess.H6x6[source.vecIndex].block<3, 3>(3 * source.blockIndices[0], 3 * source.blockIndices[1]);
                            break;
                        case LinearHessBlockType::Block9x9:
                            accBlock += barrierHess.H9x9[source.vecIndex].block<3, 3>(3 * source.blockIndices[0], 3 * source.blockIndices[1]);
                            break;
                        case LinearHessBlockType::Block12x12:
                            accBlock += barrierHess.H12x12[source.vecIndex].block<3, 3>(3 * source.blockIndices[0], 3 * source.blockIndices[1]);
                            break;
                        default:
                            BSIPC_ERROR("Unknown linear block type");
                        }
                    }
                    block.second.first = accBlock;
                }
            );
        }
        BSIPC_PROFILE_END("Accumulate Linear Hess Blocks");

        BSIPC_PROFILE_START("Create/Flatten Hash");
        // [dof] prefix describes properties of the final system.
        //       This could include both B-spline control points, and animated mesh vertices.
        using BSHessBlockCache = std::pair<Mat<3, 3>, tbb::concurrent_vector<BSHessBlockSource>>;
        tbb::concurrent_unordered_map<uint64_t, BSHessBlockCache> dofBlockSourceMap;
        BSIPC_PROFILE_END("Create/Flatten Hash");

        BSIPC_PROFILE_START("Record Active DoF Hess Block Sources");
        {
            UInt bsTrigVertEndIndex = 0;
            for (const BSTargetInfo& info : this->targets)
                bsTrigVertEndIndex += info.trigVertCnt.value();

            UInt bsDofEndIndex = this->bsVertCntSum;

            tbb::parallel_for_each(
                linearBlockSourceMap.begin(), linearBlockSourceMap.end(),
                [&](auto& block)
                {
                    std::pair<UInt, UInt> rowCol = DecompressHash(block.first);                     // Indices in the linear triangle mesh
                    UInt row = rowCol.first, col = rowCol.second;

                    // avoid multiple heap allocation, use fixed-size array instead. [row/col]DofIndicesCnt indicates the actual number of entries in the array.
                    std::array<UInt, 9> rowDofIndices, colDofIndices;
                    std::array<Float, 9> rowDofCoeffs, colDofCoeffs;
                    UInt rowDofIndicesCnt = 0, colDofIndicesCnt = 0;

                    // Propagate indices from linear triangle mesh to system DoFs
                    if (row < bsTrigVertEndIndex)                                                   // Belongs to BS surface
                    {
                        const std::vector<SourceInfo>& sources = this->linearVertexSources[row];
                        //BSIPC_ASSERT(sources.size() == 9, "Check: In [InitContactMesh] it is required to have sources.size() == 9");

                        for (UInt i = 0; i != 9; ++i)
                        {
                            //BSIPC_ASSERT(sources[i].coeff >= 0, "Coefficient associating linear triangle vertex to B-spline control point must be positive");
                            if (sources[i].coeff > 1e-5)
                            {
                                rowDofIndices[rowDofIndicesCnt] = sources[i].index;
                                rowDofCoeffs[rowDofIndicesCnt] = sources[i].coeff;
                                //BSIPC_ASSERT(rowDofIndices[rowDofIndicesCnt] < size, "Block index out of range");

                                ++rowDofIndicesCnt;
                            }
                        }
                    }
                    else if (this->includeAnimatedMeshInSystem)                     // only add non-BS DoFs when animated mesh should be included in the system
                    {
                        rowDofIndices[0] = (row - bsTrigVertEndIndex) + bsDofEndIndex;
                        rowDofCoeffs[0] = 1.;
                        //BSIPC_ASSERT(rowDofIndices[rowDofIndicesCnt] < size, "Block index out of range");

                        rowDofIndicesCnt = 1;
                    }

                    if (col < bsTrigVertEndIndex)                                                   // Belongs to BS surface
                    {
                        const std::vector<SourceInfo>& sources = this->linearVertexSources[col];
                        //BSIPC_ASSERT(sources.size() == 9, "Check: In [InitContactMesh] it is required to have sources.size() == 9");

                        for (UInt i = 0; i != 9; ++i)
                        {
                            //BSIPC_ASSERT(sources[i].coeff >= 0, "Coefficient associating linear triangle vertex to B-spline control point must be positive");
                            if (sources[i].coeff > 1e-5)
                            {
                                colDofIndices[colDofIndicesCnt] = sources[i].index;
                                colDofCoeffs[colDofIndicesCnt] = sources[i].coeff;
                                //BSIPC_ASSERT(colDofIndices[colDofIndicesCnt] < size, "Block index out of range");

                                ++colDofIndicesCnt;
                            }
                        }
                    }
                    else if (this->includeAnimatedMeshInSystem)                     // only add non-BS DoFs when animated mesh should be included in the system
                    {
                        colDofIndices[0] = (col - bsTrigVertEndIndex) + bsDofEndIndex;
                        colDofCoeffs[0] = 1.;
                        //BSIPC_ASSERT(colDofIndices[colDofIndicesCnt] < size, "Block index out of range");

                        colDofIndicesCnt = 1;
                    }

                    // Active indices pairs are rowDofIndices x colDofIndices.
                    for (UInt i = 0; i != rowDofIndicesCnt; ++i)
                        for (UInt j = 0; j != colDofIndicesCnt; ++j)
                        {
                            UInt blockRow = rowDofIndices[i], blockCol = colDofIndices[j];

                            auto insertResult = dofBlockSourceMap.insert({
                                CompressKey(blockRow, blockCol),
                                std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<BSHessBlockSource>{})
                            });

                            BSHessBlockSource source{
                                .sourceBlock = &block.second.first,
                                .coeff = rowDofCoeffs[i] * colDofCoeffs[j]
                            };
                            insertResult.first->second.second.emplace_back(source);
                        }
                }
            );
        }
        BSIPC_PROFILE_END("Record Active DoF Hess Block Sources");

        BSIPC_PEEK(dofBlockSourceMap.size());

        BSIPC_PROFILE_START("Accumulate DoF (BS + Animation Sequence) Hess Blocks");
        {
            tbb::parallel_for_each(
                dofBlockSourceMap.begin(), dofBlockSourceMap.end(),
                [&](auto& block)
                {
                    Mat<3, 3> accBlock = Mat<3, 3>::Zero();
                    for (const BSHessBlockSource& source : block.second.second)
                        accBlock += source.coeff * (*source.sourceBlock);
                    block.second.first = accBlock;
                }
            );
        }
        BSIPC_PROFILE_END("Accumulate DoF (BS + Animation Sequence) Hess Blocks");

        BSIPC_PROFILE_START("Manual Sparse Matrix Construction");
        {
            /// Pass 1: Sort all the block indices over column indices

            // outer index: col index; inner index: row index
            std::vector<std::vector<UInt>> blockRowColIndices;
            blockRowColIndices.resize(size);

            //tbb::parallel_for_each(
                //dofBlockSourceMap.begin(), dofBlockSourceMap.end(),
                //[&](auto& block)
            for (const auto& block : dofBlockSourceMap)
                {
                    std::pair<UInt, UInt> rowCol = DecompressHash(block.first);
                    UInt row = rowCol.first, col = rowCol.second;
                    //BSIPC_ASSERT(row < size && col < size, "Block index out of range");
                    blockRowColIndices[col].emplace_back(row);
                }
            //);

            /// Pass 2: Sort within each column, and record the outer starts
            tbb::parallel_for(
                0, static_cast<Int>(blockRowColIndices.size()), 1, [&](Int col)
                {
                    // use std::sort as matrix is sparse, and each (block) column has few entries
                    std::sort(blockRowColIndices[col].begin(), blockRowColIndices[col].end());
                }
            );

            hess.resize(3 * size, 3 * size);
            Int* hessOuterStarts = hess.outerIndexPtr();
            hessOuterStarts[0] = 0;
            for (UInt col = 0; col != size; ++col)
            {
                UInt nnzInThisCol = 3 * static_cast<UInt>(blockRowColIndices[col].size());
                for (UInt i = 0; i != 3; ++i)
                {
                    UInt index = 3 * col + i + 1;
                    //BSIPC_ASSERT(index > 0 && index < 3 * size + 1, "attempting to access memory out of bound");
                    hessOuterStarts[index] = hessOuterStarts[index - 1] + nnzInThisCol;
                }
            }

            /// Pass 3: Fill in the entries
            UInt nnzCount = dofBlockSourceMap.size() * 9;           // each block is 3-by-3
            //BSIPC_ASSERT(hessOuterStarts[3 * size] == nnzCount, "the last entry in [outerStarts] vec should be #nnz");
            hess.reserve(nnzCount);
            Int* hessInnerIndices = hess.innerIndexPtr();
            Float* hessValues = hess.valuePtr();

            tbb::parallel_for(
                0, static_cast<Int>(size), 1, [&](Int blockCol)
            //for (UInt blockCol = 0; blockCol != size; ++blockCol)
                {
                    UInt blockRowCnt = static_cast<UInt>(blockRowColIndices[blockCol].size());
                    for (UInt blockRowIdx = 0; blockRowIdx != blockRowCnt; ++blockRowIdx)
                    {
                        UInt blockRow = blockRowColIndices[blockCol][blockRowIdx];
                        Mat<3, 3> block = dofBlockSourceMap[CompressKey(blockRow, blockCol)].first;
                        for (UInt i = 0; i != 3; ++i)               // for each column of the block
                            for (UInt j = 0; j != 3; ++j)           // entry at (j, i)
                            {
                                UInt nnzIndex = hessOuterStarts[3 * blockCol + i] + 3 * blockRowIdx + j;
                                hessValues[nnzIndex] = block(j, i);
                                hessInnerIndices[nnzIndex] = 3 * blockRow + j;
                            }
                    }
                }
            );

            hess.resizeNonZeros(nnzCount);
            hess.makeCompressed();
        }
        BSIPC_PROFILE_END("Manual Sparse Matrix Construction");

        BSIPC_PROFILE_END("[Func] ContactBarrierHessReordered");
    }

    void Solver::ContactBarrierHessReorderedOpt(SpMat& hess)
    {
        BSIPC_PROFILE_START("[Func] ContactBarrierHessReordered");

        BSIPC_PROFILE_START("Init");
        const IPC::BHessian& barrierHess = this->contactMeshBarrierHess;
        auto CompressKey = [](UInt i, UInt j) -> uint64_t {
            return (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
            };
        auto DecompressKey = [](uint64_t key) -> std::pair<UInt, UInt> {
            UInt i = static_cast<UInt>(key >> 32);
            UInt j = static_cast<UInt>(key & 0xFFFFFFFF);
            return { i, j };
            };

        UInt size = this->bsVertCntSum;
        if (this->includeAnimatedMeshInSystem)
            size += this->animatedMeshCache.value().vertCnt;

        using LinearHessBlockCache = std::pair<Mat<3, 3>, tbb::concurrent_vector<LinearHessBlockSource>>;
        tbb::concurrent_unordered_map<uint64_t, LinearHessBlockCache> linearBlockSourceMap;

        BSIPC_PROFILE_END("Init");

        BSIPC_PROFILE_START("Skim through all vertices");
        Vec3 posLowerBound = Vec3::Constant(std::numeric_limits<Float>::max()), posUpperBound = Vec3::Constant(std::numeric_limits<Float>::min());
        const std::vector<Vec3>& contactMeshVertices = this->contactMesh.vertices;
        tbb::parallel_for(
            0, static_cast<Int>(contactMeshVertices.size()), 1, [&](UInt idx)
            {
                const Vec3& pos = contactMeshVertices[idx];
                for (UInt i = 0; i != 3; ++i)
                {
                    posLowerBound[i] = std::min(posLowerBound[i], pos[i]);
                    posUpperBound[i] = std::max(posUpperBound[i], pos[i]);
                }
            }
        );
        Vec3 gridSize = (posUpperBound - posLowerBound) / 4.;
        BSIPC_PROFILE_END("Skim through all vertices");

        BSIPC_INFO(
            "#D1Index: {}, #D2Index: {}, #D3Index: {}, #D4Index: {}",
            barrierHess.D1Index.size(), barrierHess.D2Index.size(), barrierHess.D3Index.size(), barrierHess.D4Index.size()
        );

        BSIPC_PROFILE_START("Record Linear Hess Block Sources");
        {
            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.H3x3.size()), 1, [&](UInt idx)
                {
                    uint64_t key = CompressKey(barrierHess.D1Index[idx], barrierHess.D1Index[idx]);
                    LinearHessBlockSource source{
                        .blockType = LinearHessBlockType::Block3x3,
                        .vecIndex = idx,
                        .blockIndices = {0, 0}
                    };

                    auto insertResult = linearBlockSourceMap.insert({
                        key,
                        std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{})
                        });

                    insertResult.first->second.second.emplace_back(source);
                }
            );

            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.D2Index.size()), 1, [&](UInt idx)
                {
                    for (UInt i = 0; i != 2; ++i)
                        for (UInt j = 0; j != 2; ++j)
                        {
                            uint64_t key = CompressKey(barrierHess.D2Index[idx][i], barrierHess.D2Index[idx][j]);
                            LinearHessBlockSource source{
                                .blockType = LinearHessBlockType::Block6x6,
                                .vecIndex = idx,
                                .blockIndices = {i, j}
                            };

                            auto insertResult = linearBlockSourceMap.insert({
                                key,
                                std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{})
                                });

                            insertResult.first->second.second.emplace_back(source);
                        }
                }
            );

            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.D3Index.size()), 1, [&](UInt idx)
                {
                    for (UInt i = 0; i != 3; ++i)
                        for (UInt j = 0; j != 3; ++j)
                        {
                            uint64_t key = CompressKey(barrierHess.D3Index[idx][i], barrierHess.D3Index[idx][j]);
                            LinearHessBlockSource source{
                                .blockType = LinearHessBlockType::Block9x9,
                                .vecIndex = idx,
                                .blockIndices = {i, j}
                            };

                            auto insertResult = linearBlockSourceMap.insert({
                                key,
                                std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{})
                                });

                            insertResult.first->second.second.emplace_back(source);
                        }
                }
            );

            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.D4Index.size()), 1, [&](UInt idx)
                {
                    for (UInt i = 0; i != 4; ++i)
                        for (UInt j = 0; j != 4; ++j)
                        {
                            uint64_t key = CompressKey(barrierHess.D4Index[idx][i], barrierHess.D4Index[idx][j]);
                            LinearHessBlockSource source{
                                .blockType = LinearHessBlockType::Block12x12,
                                .vecIndex = idx,
                                .blockIndices = {i, j}
                            };

                            auto insertResult = linearBlockSourceMap.insert({
                                key,
                                std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{})
                                });

                            insertResult.first->second.second.emplace_back(source);
                        }
                }
            );
        }
        BSIPC_PROFILE_END("Record Linear Hess Block Sources");

        BSIPC_INFO("Linear Block Source Map Size: {}", linearBlockSourceMap.size());

        BSIPC_PROFILE_START("Accumulate Linear Hess Blocks");
        {
            tbb::parallel_for_each(
                linearBlockSourceMap.begin(), linearBlockSourceMap.end(),
                [&](auto& block)
                {
                    Mat<3, 3> accBlock = Mat<3, 3>::Zero();
                    for (const LinearHessBlockSource& source : block.second.second)
                    {
                        switch (source.blockType)
                        {
                        case LinearHessBlockType::Block3x3:
                            accBlock += barrierHess.H3x3[source.vecIndex];
                            break;
                        case LinearHessBlockType::Block6x6:
                            accBlock += barrierHess.H6x6[source.vecIndex].block<3, 3>(3 * source.blockIndices[0], 3 * source.blockIndices[1]);
                            break;
                        case LinearHessBlockType::Block9x9:
                            accBlock += barrierHess.H9x9[source.vecIndex].block<3, 3>(3 * source.blockIndices[0], 3 * source.blockIndices[1]);
                            break;
                        case LinearHessBlockType::Block12x12:
                            accBlock += barrierHess.H12x12[source.vecIndex].block<3, 3>(3 * source.blockIndices[0], 3 * source.blockIndices[1]);
                            break;
                        default:
                            BSIPC_ERROR("Unknown linear block type");
                        }
                    }
                    block.second.first = accBlock;
                }
            );
        }
        BSIPC_PROFILE_END("Accumulate Linear Hess Blocks");

        BSIPC_PROFILE_START("Create/Flatten Hash");
        // [dof] prefix describes properties of the final system.
        //       This could include both B-spline control points, and animated mesh vertices.
        using BSHessBlockCache = std::pair<Mat<3, 3>, tbb::concurrent_vector<BSHessBlockSource>>;
        tbb::concurrent_unordered_map<uint64_t, UInt> dofBlockSourceMap;
        tbb::concurrent_unordered_map<uint64_t, tbb::spin_mutex> dofBlockMutexMap;
        std::atomic<UInt> dofBlockCnt = 0;
        std::vector<std::tuple<uint64_t, Mat<3, 3>, BSHessBlockCache>> dofBlockSourceVec;
        dofBlockSourceMap.reserve(18 * linearBlockSourceMap.size());
        dofBlockMutexMap.reserve(18 * linearBlockSourceMap.size());
        dofBlockSourceVec.resize(18 * linearBlockSourceMap.size());
        BSIPC_PROFILE_END("Create/Flatten Hash");

        BSIPC_PROFILE_START("Categorize vertices into buckets");
        std::vector<std::vector<std::pair<uint64_t, const Mat<3, 3>*>>> spatialBuckets;
        UInt gridCnt = 8;
        spatialBuckets.resize(gridCnt * gridCnt * gridCnt);
        {
            for (const auto& block : linearBlockSourceMap)
            {
                IVec3 bucketVecIdx = IVec3::Zero();
                auto [blockColIdx, blockRowIdx] = DecompressKey(block.first);

                const Vec3& pos = contactMeshVertices[blockColIdx];
                for (UInt i = 0; i != 3; ++i)
                    bucketVecIdx[i] = std::max(static_cast<Int>(std::floor((pos[i] - posLowerBound[i]) / gridSize[i])), 0);

                UInt blockFlattenedIdx = bucketVecIdx[0] + gridCnt * bucketVecIdx[1] + gridCnt * gridCnt * bucketVecIdx[2];
                //BSIPC_ASSERT(blockFlattenedIdx < 64 && blockFlattenedIdx >= 0, "invalid flattened bucket idx");
                std::pair<uint64_t, const Mat<3, 3>*> blockInfo = std::make_pair(block.first, &block.second.first);
                spatialBuckets[blockFlattenedIdx].emplace_back(blockInfo);
            }
        }
        BSIPC_PROFILE_END("Categorize vertices into buckets");

        BSIPC_PROFILE_START("Record Active DoF Hess Block Sources **");
        {
            UInt bsTrigVertEndIndex = 0;
            for (const BSTargetInfo& info : this->targets)
                bsTrigVertEndIndex += info.trigVertCnt.value();

            UInt bsDofEndIndex = this->bsVertCntSum;

            tbb::parallel_for(
                0, static_cast<Int>(spatialBuckets.size()), 1, [&](Int bucket)
                {
                    const std::vector<std::pair<uint64_t, const Mat<3, 3>*>>& bucketData = spatialBuckets[bucket];

                    tbb::parallel_for_each(
                        bucketData.begin(), bucketData.end(),
                        [&](const std::pair<uint64_t, const Mat<3, 3>*>& block)
                        //for (const auto& block : bucketData)
                        {
                            auto [row, col] = DecompressKey(block.first);                     // Indices in the linear triangle mesh

                            // avoid multiple heap allocation, use fixed-size array instead. [row/col]DofIndicesCnt indicates the actual number of entries in the array.
                            std::array<UInt, 9> rowDofIndices, colDofIndices;
                            std::array<Float, 9> rowDofCoeffs, colDofCoeffs;
                            UInt rowDofIndicesCnt = 0, colDofIndicesCnt = 0;

                            // Propagate indices from linear triangle mesh to system DoFs
                            if (row < bsTrigVertEndIndex)                                                   // Belongs to BS surface
                            {
                                const std::vector<SourceInfo>& sources = this->linearVertexSources[row];
                                //BSIPC_ASSERT(sources.size() == 9, "Check: In [InitContactMesh] it is required to have sources.size() == 9");

                                for (UInt i = 0; i != 9; ++i)
                                {
                                    //BSIPC_ASSERT(sources[i].coeff >= 0, "Coefficient associating linear triangle vertex to B-spline control point must be positive");
                                    if (sources[i].coeff > 1e-5)
                                    {
                                        rowDofIndices[rowDofIndicesCnt] = sources[i].index;
                                        rowDofCoeffs[rowDofIndicesCnt] = sources[i].coeff;
                                        //BSIPC_ASSERT(rowDofIndices[rowDofIndicesCnt] < size, "Block index out of range");

                                        ++rowDofIndicesCnt;
                                    }
                                }
                            }
                            else if (this->includeAnimatedMeshInSystem)                     // only add non-BS DoFs when animated mesh should be included in the system
                            {
                                rowDofIndices[0] = (row - bsTrigVertEndIndex) + bsDofEndIndex;
                                rowDofCoeffs[0] = 1.;
                                //BSIPC_ASSERT(rowDofIndices[rowDofIndicesCnt] < size, "Block index out of range");

                                rowDofIndicesCnt = 1;
                            }

                            if (col < bsTrigVertEndIndex)                                                   // Belongs to BS surface
                            {
                                const std::vector<SourceInfo>& sources = this->linearVertexSources[col];
                                //BSIPC_ASSERT(sources.size() == 9, "Check: In [InitContactMesh] it is required to have sources.size() == 9");

                                for (UInt i = 0; i != 9; ++i)
                                {
                                    //BSIPC_ASSERT(sources[i].coeff >= 0, "Coefficient associating linear triangle vertex to B-spline control point must be positive");
                                    if (sources[i].coeff > 1e-5)
                                    {
                                        colDofIndices[colDofIndicesCnt] = sources[i].index;
                                        colDofCoeffs[colDofIndicesCnt] = sources[i].coeff;
                                        //BSIPC_ASSERT(colDofIndices[colDofIndicesCnt] < size, "Block index out of range");

                                        ++colDofIndicesCnt;
                                    }
                                }
                            }
                            else if (this->includeAnimatedMeshInSystem)                     // only add non-BS DoFs when animated mesh should be included in the system
                            {
                                colDofIndices[0] = (col - bsTrigVertEndIndex) + bsDofEndIndex;
                                colDofCoeffs[0] = 1.;
                                //BSIPC_ASSERT(colDofIndices[colDofIndicesCnt] < size, "Block index out of range");

                                colDofIndicesCnt = 1;
                            }

                            // Active indices pairs are rowDofIndices x colDofIndices.
                            for (UInt i = 0; i != rowDofIndicesCnt; ++i)
                                for (UInt j = 0; j != colDofIndicesCnt; ++j)
                                {
                                    UInt blockRow = rowDofIndices[i], blockCol = colDofIndices[j];
                                    dofBlockMutexMap[CompressKey(blockRow, blockCol)];      // ensure mutex exists

                                    auto insertResult = dofBlockSourceMap.insert({
                                        CompressKey(blockRow, blockCol),
                                        std::numeric_limits<UInt>::max()
                                        });

                                    UInt& vecIdx = insertResult.first->second;

                                    dofBlockMutexMap[CompressKey(blockRow, blockCol)].lock();
                                    if (vecIdx == -1)
                                        vecIdx = dofBlockCnt.fetch_add(1);
                                    dofBlockMutexMap[CompressKey(blockRow, blockCol)].unlock();

                                    BSHessBlockSource source{
                                        .sourceBlock = block.second,
                                        .coeff = rowDofCoeffs[i] * colDofCoeffs[j]
                                    };

                                    std::get<2>(dofBlockSourceVec[vecIdx]).second.emplace_back(source);
                                    std::get<0>(dofBlockSourceVec[vecIdx]) = CompressKey(blockRow, blockCol);
                                }
                        }
                    );
                }
            );
        }
        BSIPC_PROFILE_END("Record Active DoF Hess Block Sources **");

        BSIPC_INFO("dof Block Source Map Size: {}", dofBlockSourceMap.size());

        BSIPC_PROFILE_START("Accumulate DoF (BS + Animation Sequence) Hess Blocks");
        {
            // Pre-extract to vector for better cache locality
            //dofBlockSourceVec.reserve(dofBlockSourceMap.size());
            //for (auto& block : dofBlockSourceMap)
            //    dofBlockSourceVec.emplace_back(block.first, Mat<3, 3>::Zero(), &block.second);

            tbb::parallel_for(
                0, static_cast<Int>(dofBlockCnt.load()), 1, 
                [&](Int i)
                {
                    auto& [key, accBlock, cache] = dofBlockSourceVec[i];
                    for (const BSHessBlockSource& source : cache.second)
                        accBlock += source.coeff * (*source.sourceBlock);
                }
            );

            //tbb::parallel_for_each(
            //    dofBlockSourceMap.begin(), dofBlockSourceMap.end(),
            //    [&](auto& block)
            //    {
            //        Mat<3, 3> accBlock = Mat<3, 3>::Zero();
            //        for (const BSHessBlockSource& source : block.second.second)
            //            accBlock += source.coeff * (*source.sourceBlock);
            //        block.second.first = accBlock;
            //    }
            //);
        }
        BSIPC_PROFILE_END("Accumulate DoF (BS + Animation Sequence) Hess Blocks");

        BSIPC_INFO("dofBlockCnt: {}, #ofMap: {}", dofBlockCnt.load(), dofBlockSourceMap.size());

        BSIPC_PROFILE_START("Manual Sparse Matrix Construction");
        {
            /// Pass 1: Categorize all the block indices over column indices

            // outer index: col index; inner index: row index
            std::vector<std::vector<std::pair<UInt, Mat<3, 3>*>>> blockRowColIndices;
            blockRowColIndices.resize(size);

            //tbb::parallel_for_each(
                //dofBlockSourceMap.begin(), dofBlockSourceMap.end(),
                //[&](auto& block)
            //for (const auto& block : dofBlockSourceMap)
            //{
            //    std::pair<UInt, UInt> rowCol = DecompressKey(block.first);
            //    UInt row = rowCol.first, col = rowCol.second;
            //    //BSIPC_ASSERT(row < size && col < size, "Block index out of range");
            //    blockRowColIndices[col].emplace_back(row);
            //}
            //);

            //tbb::parallel_for(
                //0, static_cast<Int>(dofBlockSourceVec.size()), 1,
                //[&](Int i)
            for (UInt i = 0; i != dofBlockCnt.load(); ++i)
                {
                    auto& [key, blockHess, _sources] = dofBlockSourceVec[i];
                    std::pair<UInt, UInt> rowCol = DecompressKey(key);
                    UInt row = rowCol.first, col = rowCol.second;
                    //BSIPC_ASSERT(row < size && col < size, "Block index out of range");
                    blockRowColIndices[col].emplace_back(row, &blockHess);
                }
            //);

            /// Pass 2: Sort within each column, and record the outer starts
            tbb::parallel_for(
                0, static_cast<Int>(blockRowColIndices.size()), 1, 
                [&](Int col)
                {
                    // use std::sort instead of tbb::parallel_sort as matrix is sparse, and each (block) column has few entries
                    std::sort(
                        blockRowColIndices[col].begin(), blockRowColIndices[col].end(), 
                        [](const std::pair<UInt, Mat<3, 3>*>& a, const std::pair<UInt, Mat<3, 3>*>& b) { return a.first < b.first; }
                    );
                }
            );

            UInt hessSize = 3 * size;
            hess.resize(hessSize, hessSize);
            UInt nnzCount = dofBlockSourceMap.size() * 9;           // each block is 3-by-3
            hess.reserve(nnzCount);
            Int* hessOuterStarts = hess.outerIndexPtr();
            hessOuterStarts[0] = 0;
            for (UInt col = 0; col != size; ++col)
            {
                UInt nnzInThisCol = 3 * static_cast<UInt>(blockRowColIndices[col].size());
                for (UInt i = 0; i != 3; ++i)
                {
                    UInt index = 3 * col + i + 1;
                    //BSIPC_ASSERT(index > 0 && index < 3 * size + 1, "attempting to access memory out of bound");
                    hessOuterStarts[index] = hessOuterStarts[index - 1] + nnzInThisCol;
                }
            }
            BSIPC_ASSERT(hessOuterStarts[3 * size] == nnzCount, "the last entry in [outerStarts] vec should be #nnz");

            /// Pass 3: Fill in the entries
            Int* hessInnerIndices = hess.innerIndexPtr();
            Float* hessValues = hess.valuePtr();

            tbb::parallel_for(
                0, static_cast<Int>(size), 1, [&](Int blockCol)
                //for (UInt blockCol = 0; blockCol != size; ++blockCol)
                {
                    UInt blockRowCnt = static_cast<UInt>(blockRowColIndices[blockCol].size());
                    for (UInt blockRowIdx = 0; blockRowIdx != blockRowCnt; ++blockRowIdx)
                    {
                        auto& [blockRow, block] = blockRowColIndices[blockCol][blockRowIdx];
                        for (UInt i = 0; i != 3; ++i)               // for each column of the block
                            for (UInt j = 0; j != 3; ++j)           // entry at (j, i)
                            {
                                UInt nnzIndex = hessOuterStarts[3 * blockCol + i] + 3 * blockRowIdx + j;
                                hessValues[nnzIndex] = (*block)(j, i);
                                hessInnerIndices[nnzIndex] = 3 * blockRow + j;
                            }
                    }
                }
            );

            hess.resizeNonZeros(nnzCount);
            hess.makeCompressed();
        }
        BSIPC_PROFILE_END("Manual Sparse Matrix Construction");

        BSIPC_PROFILE_END("[Func] ContactBarrierHessReordered");
    }

    void Solver::ContactBarrierHessReordered(SpMat& hess)
    {
        BSIPC_PROFILE_START("[Func] ContactBarrierHessReordered");

        BSIPC_PROFILE_START("Init");
        const IPC::BHessian& barrierHess = this->contactMeshBarrierHess;
        auto CompressKey = [](UInt i, UInt j) -> uint64_t {
            return (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
            };
        auto DecompressKey = [](uint64_t key) -> std::pair<UInt, UInt> {
            UInt i = static_cast<UInt>(key >> 32);
            UInt j = static_cast<UInt>(key & 0xFFFFFFFF);
            return { i, j };
            };

        UInt size = this->bsVertCntSum;
        if (this->includeAnimatedMeshInSystem)
            size += this->animatedMeshCache.value().vertCnt;

        using LinearHessBlockCache = std::pair<Mat<3, 3>, tbb::concurrent_vector<LinearHessBlockSource>>;
        tbb::concurrent_unordered_map<uint64_t, LinearHessBlockCache> linearBlockSourceMap;

        BSIPC_PROFILE_END("Init");

        BSIPC_PROFILE_START("Skim through all vertices");
        Vec3 posLowerBound = Vec3::Constant(std::numeric_limits<Float>::max()), posUpperBound = Vec3::Constant(std::numeric_limits<Float>::min());
        const std::vector<Vec3>& contactMeshVertices = this->contactMesh.vertices;
        tbb::parallel_for(
            0, static_cast<Int>(contactMeshVertices.size()), 1, [&](UInt idx)
            {
                const Vec3& pos = contactMeshVertices[idx];
                for (UInt i = 0; i != 3; ++i)
                {
                    posLowerBound[i] = std::min(posLowerBound[i], pos[i]);
                    posUpperBound[i] = std::max(posUpperBound[i], pos[i]);
                }
            }
        );
        Vec3 gridSize = (posUpperBound - posLowerBound) / 4.;
        BSIPC_PROFILE_END("Skim through all vertices");

        BSIPC_INFO(
            "#D1Index: {}, #D2Index: {}, #D3Index: {}, #D4Index: {}",
            barrierHess.D1Index.size(), barrierHess.D2Index.size(), barrierHess.D3Index.size(), barrierHess.D4Index.size()
        );

        BSIPC_PROFILE_START("Record Linear Hess Block Sources");
        {
            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.H3x3.size()), 1, [&](UInt idx)
                {
                    uint64_t key = CompressKey(barrierHess.D1Index[idx], barrierHess.D1Index[idx]);
                    LinearHessBlockSource source{
                        .blockType = LinearHessBlockType::Block3x3,
                        .vecIndex = idx,
                        .blockIndices = {0, 0}
                    };

                    auto insertResult = linearBlockSourceMap.insert({
                        key,
                        std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{})
                        });

                    insertResult.first->second.second.emplace_back(source);
                }
            );

            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.D2Index.size()), 1, [&](UInt idx)
                {
                    for (UInt i = 0; i != 2; ++i)
                        for (UInt j = 0; j != 2; ++j)
                        {
                            uint64_t key = CompressKey(barrierHess.D2Index[idx][i], barrierHess.D2Index[idx][j]);
                            LinearHessBlockSource source{
                                .blockType = LinearHessBlockType::Block6x6,
                                .vecIndex = idx,
                                .blockIndices = {i, j}
                            };

                            auto insertResult = linearBlockSourceMap.insert({
                                key,
                                std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{})
                                });

                            insertResult.first->second.second.emplace_back(source);
                        }
                }
            );

            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.D3Index.size()), 1, [&](UInt idx)
                {
                    for (UInt i = 0; i != 3; ++i)
                        for (UInt j = 0; j != 3; ++j)
                        {
                            uint64_t key = CompressKey(barrierHess.D3Index[idx][i], barrierHess.D3Index[idx][j]);
                            LinearHessBlockSource source{
                                .blockType = LinearHessBlockType::Block9x9,
                                .vecIndex = idx,
                                .blockIndices = {i, j}
                            };

                            auto insertResult = linearBlockSourceMap.insert({
                                key,
                                std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{})
                                });

                            insertResult.first->second.second.emplace_back(source);
                        }
                }
            );

            tbb::parallel_for(
                0, static_cast<Int>(barrierHess.D4Index.size()), 1, [&](UInt idx)
                {
                    for (UInt i = 0; i != 4; ++i)
                        for (UInt j = 0; j != 4; ++j)
                        {
                            uint64_t key = CompressKey(barrierHess.D4Index[idx][i], barrierHess.D4Index[idx][j]);
                            LinearHessBlockSource source{
                                .blockType = LinearHessBlockType::Block12x12,
                                .vecIndex = idx,
                                .blockIndices = {i, j}
                            };

                            auto insertResult = linearBlockSourceMap.insert({
                                key,
                                std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<LinearHessBlockSource>{})
                                });

                            insertResult.first->second.second.emplace_back(source);
                        }
                }
            );
        }
        BSIPC_PROFILE_END("Record Linear Hess Block Sources");

        BSIPC_INFO("Linear Block Source Map Size: {}", linearBlockSourceMap.size());

        BSIPC_PROFILE_START("Accumulate Linear Hess Blocks");
        {
            tbb::parallel_for_each(
                linearBlockSourceMap.begin(), linearBlockSourceMap.end(),
                [&](auto& block)
                {
                    Mat<3, 3> accBlock = Mat<3, 3>::Zero();
                    for (const LinearHessBlockSource& source : block.second.second)
                    {
                        switch (source.blockType)
                        {
                        case LinearHessBlockType::Block3x3:
                            accBlock += barrierHess.H3x3[source.vecIndex];
                            break;
                        case LinearHessBlockType::Block6x6:
                            accBlock += barrierHess.H6x6[source.vecIndex].block<3, 3>(3 * source.blockIndices[0], 3 * source.blockIndices[1]);
                            break;
                        case LinearHessBlockType::Block9x9:
                            accBlock += barrierHess.H9x9[source.vecIndex].block<3, 3>(3 * source.blockIndices[0], 3 * source.blockIndices[1]);
                            break;
                        case LinearHessBlockType::Block12x12:
                            accBlock += barrierHess.H12x12[source.vecIndex].block<3, 3>(3 * source.blockIndices[0], 3 * source.blockIndices[1]);
                            break;
                        default:
                            BSIPC_ERROR("Unknown linear block type");
                        }
                    }
                    block.second.first = accBlock;
                }
            );
        }
        BSIPC_PROFILE_END("Accumulate Linear Hess Blocks");

        BSIPC_PROFILE_START("Create/Flatten Hash");
        // [dof] prefix describes properties of the final system.
        //       This could include both B-spline control points, and animated mesh vertices.
        using BSHessBlockCache = std::pair<Mat<3, 3>, tbb::concurrent_vector<BSHessBlockSource>>;
        tbb::concurrent_unordered_map<uint64_t, BSHessBlockCache> dofBlockSourceMap;
        BSIPC_PROFILE_END("Create/Flatten Hash");

        BSIPC_PROFILE_START("Categorize vertices into buckets");
        std::vector<std::vector<std::pair<uint64_t, const Mat<3, 3>*>>> spatialBuckets;
        UInt gridCnt = 8;
        spatialBuckets.resize(gridCnt* gridCnt* gridCnt);
        {
            for (const auto& block : linearBlockSourceMap)
            {
                IVec3 bucketVecIdx = IVec3::Zero();
                auto [blockColIdx, blockRowIdx] = DecompressKey(block.first);

                const Vec3& pos = contactMeshVertices[blockColIdx];
                for (UInt i = 0; i != 3; ++i)
                    bucketVecIdx[i] = std::max(static_cast<Int>(std::floor((pos[i] - posLowerBound[i]) / gridSize[i])), 0);

                UInt blockFlattenedIdx = bucketVecIdx[0] + gridCnt * bucketVecIdx[1] + gridCnt * gridCnt * bucketVecIdx[2];
                //BSIPC_ASSERT(blockFlattenedIdx < 64 && blockFlattenedIdx >= 0, "invalid flattened bucket idx");
                std::pair<uint64_t, const Mat<3, 3>*> blockInfo = std::make_pair(block.first, &block.second.first);
                spatialBuckets[blockFlattenedIdx].emplace_back(blockInfo);
            }
        }
        BSIPC_PROFILE_END("Categorize vertices into buckets");

        BSIPC_PROFILE_START("Record Active DoF Hess Block Sources **");
        {
            UInt bsTrigVertEndIndex = 0;
            for (const BSTargetInfo& info : this->targets)
                bsTrigVertEndIndex += info.trigVertCnt.value();

            UInt bsDofEndIndex = this->bsVertCntSum;

            tbb::parallel_for(
                0, static_cast<Int>(spatialBuckets.size()), 1, [&](Int bucket)
                {
                    const std::vector<std::pair<uint64_t, const Mat<3, 3>*>>& bucketData = spatialBuckets[bucket];

                    tbb::parallel_for_each(
                        bucketData.begin(), bucketData.end(),
                        [&](const std::pair<uint64_t, const Mat<3, 3>*>& block)
                        //for (const auto& block : bucketData)
                        {
                            auto [row, col] = DecompressKey(block.first);                     // Indices in the linear triangle mesh

                            // avoid multiple heap allocation, use fixed-size array instead. [row/col]DofIndicesCnt indicates the actual number of entries in the array.
                            std::array<UInt, 9> rowDofIndices, colDofIndices;
                            std::array<Float, 9> rowDofCoeffs, colDofCoeffs;
                            UInt rowDofIndicesCnt = 0, colDofIndicesCnt = 0;

                            // Propagate indices from linear triangle mesh to system DoFs
                            if (row < bsTrigVertEndIndex)                                                   // Belongs to BS surface
                            {
                                const std::vector<SourceInfo>& sources = this->linearVertexSources[row];
                                //BSIPC_ASSERT(sources.size() == 9, "Check: In [InitContactMesh] it is required to have sources.size() == 9");

                                for (UInt i = 0; i != 9; ++i)
                                {
                                    //BSIPC_ASSERT(sources[i].coeff >= 0, "Coefficient associating linear triangle vertex to B-spline control point must be positive");
                                    if (sources[i].coeff > 1e-5)
                                    {
                                        rowDofIndices[rowDofIndicesCnt] = sources[i].index;
                                        rowDofCoeffs[rowDofIndicesCnt] = sources[i].coeff;
                                        //BSIPC_ASSERT(rowDofIndices[rowDofIndicesCnt] < size, "Block index out of range");

                                        ++rowDofIndicesCnt;
                                    }
                                }
                            }
                            else if (this->includeAnimatedMeshInSystem)                     // only add non-BS DoFs when animated mesh should be included in the system
                            {
                                rowDofIndices[0] = (row - bsTrigVertEndIndex) + bsDofEndIndex;
                                rowDofCoeffs[0] = 1.;
                                //BSIPC_ASSERT(rowDofIndices[rowDofIndicesCnt] < size, "Block index out of range");

                                rowDofIndicesCnt = 1;
                            }

                            if (col < bsTrigVertEndIndex)                                                   // Belongs to BS surface
                            {
                                const std::vector<SourceInfo>& sources = this->linearVertexSources[col];
                                //BSIPC_ASSERT(sources.size() == 9, "Check: In [InitContactMesh] it is required to have sources.size() == 9");

                                for (UInt i = 0; i != 9; ++i)
                                {
                                    //BSIPC_ASSERT(sources[i].coeff >= 0, "Coefficient associating linear triangle vertex to B-spline control point must be positive");
                                    if (sources[i].coeff > 1e-5)
                                    {
                                        colDofIndices[colDofIndicesCnt] = sources[i].index;
                                        colDofCoeffs[colDofIndicesCnt] = sources[i].coeff;
                                        //BSIPC_ASSERT(colDofIndices[colDofIndicesCnt] < size, "Block index out of range");

                                        ++colDofIndicesCnt;
                                    }
                                }
                            }
                            else if (this->includeAnimatedMeshInSystem)                     // only add non-BS DoFs when animated mesh should be included in the system
                            {
                                colDofIndices[0] = (col - bsTrigVertEndIndex) + bsDofEndIndex;
                                colDofCoeffs[0] = 1.;
                                //BSIPC_ASSERT(colDofIndices[colDofIndicesCnt] < size, "Block index out of range");

                                colDofIndicesCnt = 1;
                            }

                            // Active indices pairs are rowDofIndices x colDofIndices.
                            for (UInt i = 0; i != rowDofIndicesCnt; ++i)
                                for (UInt j = 0; j != colDofIndicesCnt; ++j)
                                {
                                    UInt blockRow = rowDofIndices[i], blockCol = colDofIndices[j];

                                    auto insertResult = dofBlockSourceMap.insert({
                                        CompressKey(blockRow, blockCol),
                                        std::make_pair(Mat<3, 3>::Zero(), tbb::concurrent_vector<BSHessBlockSource>{})
                                        });

                                    BSHessBlockSource source{
                                        .sourceBlock = block.second,
                                        .coeff = rowDofCoeffs[i] * colDofCoeffs[j]
                                    };
                                    insertResult.first->second.second.emplace_back(source);
                                }
                        }
                    );
                }
            );
        }
        BSIPC_PROFILE_END("Record Active DoF Hess Block Sources **");

        BSIPC_INFO("dof Block Source Map Size: {}", dofBlockSourceMap.size());

        BSIPC_PROFILE_START("Accumulate DoF (BS + Animation Sequence) Hess Blocks");
        {
            tbb::parallel_for_each(
                dofBlockSourceMap.begin(), dofBlockSourceMap.end(),
                [&](auto& block)
                {
                    Mat<3, 3> accBlock = Mat<3, 3>::Zero();
                    for (const BSHessBlockSource& source : block.second.second)
                        accBlock += source.coeff * (*source.sourceBlock);
                    block.second.first = accBlock;
                }
            );
        }
        BSIPC_PROFILE_END("Accumulate DoF (BS + Animation Sequence) Hess Blocks");

        BSIPC_PROFILE_START("Manual Sparse Matrix Construction");
        {
            /// Pass 1: Categorize all the block indices over column indices

            // outer index: col index; inner index: row index
            std::vector<std::vector<UInt>> blockRowColIndices;
            blockRowColIndices.resize(size);

            //tbb::parallel_for_each(
                //dofBlockSourceMap.begin(), dofBlockSourceMap.end(),
                //[&](auto& block)
            for (const auto& block : dofBlockSourceMap)
            {
                std::pair<UInt, UInt> rowCol = DecompressKey(block.first);
                UInt row = rowCol.first, col = rowCol.second;
                //BSIPC_ASSERT(row < size && col < size, "Block index out of range");
                blockRowColIndices[col].emplace_back(row);
            }
            //);

            /// Pass 2: Sort within each column, and record the outer starts
            tbb::parallel_for(
                0, static_cast<Int>(blockRowColIndices.size()), 1, [&](Int col)
                {
                    // use std::sort instead of tbb::parallel_sort as matrix is sparse, and each (block) column has few entries
                    std::sort(blockRowColIndices[col].begin(), blockRowColIndices[col].end());
                }
            );

            UInt hessSize = 3 * size;
            hess.resize(hessSize, hessSize);
            UInt nnzCount = dofBlockSourceMap.size() * 9;           // each block is 3-by-3
            hess.reserve(nnzCount);
            Int* hessOuterStarts = hess.outerIndexPtr();
            hessOuterStarts[0] = 0;
            for (UInt col = 0; col != size; ++col)
            {
                UInt nnzInThisCol = 3 * static_cast<UInt>(blockRowColIndices[col].size());
                for (UInt i = 0; i != 3; ++i)
                {
                    UInt index = 3 * col + i + 1;
                    //BSIPC_ASSERT(index > 0 && index < 3 * size + 1, "attempting to access memory out of bound");
                    hessOuterStarts[index] = hessOuterStarts[index - 1] + nnzInThisCol;
                }
            }
            BSIPC_ASSERT(hessOuterStarts[3 * size] == nnzCount, "the last entry in [outerStarts] vec should be #nnz");

            /// Pass 3: Fill in the entries
            Int* hessInnerIndices = hess.innerIndexPtr();
            Float* hessValues = hess.valuePtr();

            tbb::parallel_for(
                0, static_cast<Int>(size), 1, [&](Int blockCol)
                //for (UInt blockCol = 0; blockCol != size; ++blockCol)
                {
                    UInt blockRowCnt = static_cast<UInt>(blockRowColIndices[blockCol].size());
                    for (UInt blockRowIdx = 0; blockRowIdx != blockRowCnt; ++blockRowIdx)
                    {
                        UInt blockRow = blockRowColIndices[blockCol][blockRowIdx];
                        Mat<3, 3> block = dofBlockSourceMap[CompressKey(blockRow, blockCol)].first;
                        for (UInt i = 0; i != 3; ++i)               // for each column of the block
                            for (UInt j = 0; j != 3; ++j)           // entry at (j, i)
                            {
                                UInt nnzIndex = hessOuterStarts[3 * blockCol + i] + 3 * blockRowIdx + j;
                                hessValues[nnzIndex] = block(j, i);
                                hessInnerIndices[nnzIndex] = 3 * blockRow + j;
                            }
                    }
                }
            );

            hess.resizeNonZeros(nnzCount);
            hess.makeCompressed();
        }
        BSIPC_PROFILE_END("Manual Sparse Matrix Construction");

        BSIPC_PROFILE_END("[Func] ContactBarrierHessReordered");
    }

    void Solver::_ContactBarrierHessTriplet(SpMat& hess)
    {
        BSIPC_PROFILE_START("Hess Prep");

        // No need to take into consideration the gravity term as its Hessian is zero
        UInt size = this->bsVertCntSum;

        UInt trigMeshVertCnt = 0, quadPointsCnt = 0, patchCnt = 0;
        for (const BSTargetInfo& info : this->targets)
        {
            quadPointsCnt += info.target->GetSsQuadPoints().size();
            trigMeshVertCnt += info.trigVertCnt.value();
            patchCnt += info.target->UDomainMax() * info.target->VDomainMax();
        }

        // Filter the entries in the barrier hessian s.t. only those related to the B-spline surface are involved.
        const SpMatData barrierHessData = this->contactMeshBarrierHess.toTriplets(this->contactMesh.boundaryTypes);

        tbb::concurrent_vector<SpMatEntry> bsBarrierHessData;

        UInt hessEntryLimit = 3 * trigMeshVertCnt;
        if (includeAnimatedMeshInSystem)
            hessEntryLimit += 3 * this->animatedMeshCache.value().vertCnt;
        tbb::parallel_for(
            0, static_cast<Int>(barrierHessData.size()), 1, [&](Int i)
            {
                const Eigen::Triplet<Float>& entry = barrierHessData[i];
                // Convert from global trig mesh indices to local trig mesh indices, as if this is the case where there is only one [BSSurface].
                if (entry.row() < static_cast<Int>(hessEntryLimit) && entry.col() < static_cast<Int>(hessEntryLimit))
                    bsBarrierHessData.push_back(SpMatEntry(entry.row(), entry.col(), entry.value()));
            }
        );

        BSIPC_PROFILE_END("Hess Prep");

        BSIPC_PROFILE_START("Collect Linear Hess");

        // Assemble the Hessian on the linear triangle, and repopulate into a vector.
        // This step is supposed to condense entries with multiple elements into one.
        UInt linearHessSize = 3 * trigMeshVertCnt;
        if (includeAnimatedMeshInSystem)
            linearHessSize += 3 * this->animatedMeshCache.value().vertCnt;
        SpMat linearHess(linearHessSize, linearHessSize);
        linearHess.setFromTriplets(bsBarrierHessData.begin(), bsBarrierHessData.end());

        tbb::concurrent_vector<SpMatEntry> condensedLinearHessEntries;
        tbb::concurrent_vector<SpMatEntry> bsAnimHessEntries;                       // row index: bs; col index: animated
        tbb::concurrent_vector<SpMatEntry> animBsHessEntries;                       // row index: animated; col index: bs
        tbb::concurrent_vector<SpMatEntry> animAnimHessEntries;                     // row index: animated; col index: animated

        const Int* outer = linearHess.outerIndexPtr();      // size = cols + 1
        const Int* inner = linearHess.innerIndexPtr();      // row indices
        const Float* values = linearHess.valuePtr();        // actual values
        Int outerSize = linearHess.outerSize();             // number of columns

        if (includeAnimatedMeshInSystem)
        {
            UInt bsHessLimit = 3 * trigMeshVertCnt;
            tbb::parallel_for(
                0, outerSize, [&](Int col)
                //for (Int col = 0; col < outerSize; ++col)
                {
                    for (Int idx = outer[col]; idx < outer[col + 1]; ++idx)
                    {
                        Int row = inner[idx];
                        Float val = values[idx];
                        if (row < bsHessLimit && col < bsHessLimit)
                            condensedLinearHessEntries.emplace_back(row, col, val);
                        else if (row < bsHessLimit && col >= bsHessLimit)
                            bsAnimHessEntries.emplace_back(row, col, val);
                        else if (row >= bsHessLimit && col < bsHessLimit)
                            animBsHessEntries.emplace_back(row, col, val);
                        else
                            animAnimHessEntries.emplace_back(row, col, val);
                    }
                }
            );
        }
        else
        {
            tbb::parallel_for(
                0, outerSize, [&](Int col)
                {
                    for (Int idx = outer[col]; idx < outer[col + 1]; ++idx)
                    {
                        Int row = inner[idx];
                        Float val = values[idx];
                        condensedLinearHessEntries.emplace_back(row, col, val);
                    }
                }
            );
        }

        BSIPC_PROFILE_END("Collect Linear Hess");

        //Profiler::Check("Collect Linear Hess");
        //BSIPC_TRACE("{}: {}", "Collect Linear Hess", Profiler::GetTime("Collect Linear Hess").value());
        //Profiler::Stop("Collect Linear Hess");

        UInt bsHessEntryCnt = 81 * static_cast<UInt>(condensedLinearHessEntries.size());
        UInt bsAnimHessEntryCnt = 9 * static_cast<UInt>(bsAnimHessEntries.size() + animBsHessEntries.size());
        UInt animAnimHessEntryCnt = static_cast<UInt>(animAnimHessEntries.size());

        // Each Hessian entry in the triangular mesh can affect at most 9 * 9 = 81 entries in the global Hessian
        // If the pair is one-side (bs - animated), then only the bs-side would expand to 9 entries
        const UInt barrierHessEntryCnt = bsHessEntryCnt + bsAnimHessEntryCnt + animAnimHessEntryCnt;
        //const UInt barrierHessEntryCnt = bsHessEntryCnt + bsAnimHessEntryCnt;

        this->hessEntries.clear();
        this->hessEntries.resize(barrierHessEntryCnt);

        {
            BSIPC_PROFILE_SCOPE(
                "Barrier Hess",
                //this->FillLinearHess(linearHessView, this->hessEntries, inertiaHessEntryCnt + potHessEntryCnt);
                this->FillLinearHess(condensedLinearHessEntries, this->hessEntries, 0);

            if (includeAnimatedMeshInSystem)
                this->FillAnimatedLinearHess(
                    bsAnimHessEntries, animBsHessEntries, animAnimHessEntries,
                    this->hessEntries, bsHessEntryCnt
                );
            )
        }

        if (this->includeAnimatedMeshInSystem)
            size += this->animatedMeshCache.value().vertCnt;

        hess.resize(3 * size, 3 * size);
        BSIPC_PEEK(this->hessEntries.size());
        BSIPC_PROFILE_START("Set From Triplets");
        hess.setFromTriplets(this->hessEntries.begin(), this->hessEntries.end());
        BSIPC_PROFILE_END("Set From Triplets");
    }

#if defined BSIPC_NUMERIC_TEST
    void Solver::BdHessNumericTest(const SpMatData& bdData)
    {
        UInt size = this->bsVertCntSum;

        SpMat BdHess(3 * size, 3 * size);
        BdHess.setFromTriplets(bdData.begin(), bdData.end());
        
        DMat bdHessD = SparseToDense(BdHess);
        DMat numericalBdHessD = this->config.timestep * this->config.timestep * Energy::NumericalBdHess(this->config, this->targets);
        
        for (UInt i = 0; i != 3 * size; ++i)
            for (UInt j = 0; j != 3 * size; ++j)
                if (std::abs(numericalBdHessD(i, j)) < 1e-8)
                    numericalBdHessD(i, j) = 0;
        
        this->AppendToLog("===== Hessian Assembly Test =====\n\n");
        this->AppendToLog("Analytic Bd Hess: \n" + ToStr(bdHessD) + "\n");
        this->AppendToLog("Numerical Bd Hess: \n" + ToStr(numericalBdHessD) + "\n");
        this->AppendToLog("\n<<<<< End Test >>>>>\n");
    }

    void Solver::SsGradBlockNumericTest(const QuadPoint& quadPoint, std::vector<UInt> indices, DMat blockGrad, const BSTargetInfo& info)
    {
        DMat numericalGrad = Energy::NumericalSsGradBlockAt(quadPoint, indices, this->config, info.target);
        DMat transposed = blockGrad.transpose();
        DMat numericalTransposed = numericalGrad.transpose();
        
        this->AppendToLog("===== Stretching Gradient Test =====\n\n");
        this->AppendToLog("Analytic Stretching gradient at (" + std::to_string(quadPoint.U()) + ", " + std::to_string(quadPoint.V()) + "): \n" + 
            ToStr(transposed) + "\n");
        this->AppendToLog("Numerical Stretching gradient at (" + std::to_string(quadPoint.U()) + ", " + std::to_string(quadPoint.V()) + "): \n" + 
            ToStr(numericalTransposed) + "\n");
        this->AppendToLog("\n<<<<< End Test >>>>>\n");
    }

    void Solver::SsHessBlockNumericTest(const QuadPoint & quadPoint, std::vector<UInt> indices, DMat blockHess, const BSTargetInfo & info)
    {
        DMat numericalHess = Energy::NumericalSsHessBlockAt(quadPoint, indices, this->config, info.target);
    
        this->AppendToLog("===== Stretching Hessian Test =====\n\n");
        this->AppendToLog("Analytic Stretching Hessian at (" + std::to_string(quadPoint.U()) + ", " + std::to_string(quadPoint.V()) + "): \n" +
            ToStr(blockHess) + "\n");
        this->AppendToLog("Numerical Stretching Hessian at (" + std::to_string(quadPoint.U()) + ", " + std::to_string(quadPoint.V()) + "): \n" +
            ToStr(numericalHess) + "\n");
        this->AppendToLog("\n<<<<< End Test >>>>>\n");
    }

    void Solver::BdGradBlockNumericTest(const QuadPoint& quadPoint, std::vector<UInt> indices, DMat blockGrad, const BSTargetInfo& info)
    {
        DMat numericalGrad = Energy::NumericalBdGradBlockAt(quadPoint, indices, this->config, info.target);
        DMat transposed = blockGrad.transpose();
        DMat numericalTransposed = numericalGrad.transpose();
        DMat diff = transposed - numericalTransposed;
        Float u = quadPoint.U(), v = quadPoint.V();

        for (UInt i = 0; i != diff.cols(); ++i)
            if (std::abs(diff(0, i)) < 1e-4)
                diff(0, i) = 0;

        this->AppendToLog("===== Bending Gradient Test =====\n\n");
        this->AppendToLog("Analytic Bending gradient at (" + std::to_string(u) + ", " + std::to_string(v) + "): \n" + ToStr(transposed) + "\n");
        this->AppendToLog("Numerical Bending gradient at (" + std::to_string(u) + ", " + std::to_string(v) + "): \n" + ToStr(numericalTransposed) + "\n");
        this->AppendToLog("Diff: " + ToStr(diff) + "\n");
        this->AppendToLog("\n<<<<< End Test >>>>>\n");
    }

    void Solver::BdHessBlockNumericTest(const QuadPoint& quadPoint, std::vector<UInt> indices, DMat blockHess, const BSTargetInfo& info)
    {
        // WARNING need to remove SPDProjection for bending hessian before numeric test
        DMat numericalHess = Energy::NumericalBdHessBlockAt(quadPoint, indices, this->config, info.target);
        Float u = quadPoint.U(), v = quadPoint.V();
        
        this->AppendToLog("===== Bending Hessian Test =====\n\n");
        this->AppendToLog("Analytic Bending hessian at (" + std::to_string(u) + ", " + std::to_string(v) + "): \n" + ToStr(blockHess) + "\n");
        this->AppendToLog("Numerical Bending hessian at (" + std::to_string(u) + ", " + std::to_string(v) + "): \n" + ToStr(numericalHess) + "\n");
        this->AppendToLog("\n<<<<< End Test >>>>>\n");
    }
#endif

#if defined BSIPC_NUMERIC_TEST
    void Solver::BarrierGradHessNumericTest()
    {
        IPC::BuildCollisionSets(this->contactMesh, this->spatialHash, this->ground, true);
        Float barrierVal = IPC::ComputeBarrierEnergyVal(this->contactMesh, this->ground, this->contactStiffness);

        BSIPC_PEEK(barrierVal);

        if (barrierVal > 1e-8)
        {
            UInt trigSize = static_cast<UInt>(this->trigVertCntSum);

            //DMat numericBarrierGrad_Trig = this->NumericalBarrierGrad_Trig();
            //DMat numericBarrierHess_Trig = this->NumericalBarrierHess_Trig();

            //{
            //    this->AppendToLog("===== Barrier Gradient Test [Trig] =====\n\n");
            //    DMat barGradTrigT = numericBarrierGrad_Trig.transpose();
            //    this->AppendToLog("Numerical Barrier Grad [Trig]: \n" + ToStr(barGradTrigT) + "\n");

            //    this->AppendToLog("Analytic Barrier Grad [Trig]: \n");
            //    DMat analyticBarrierGrad = DMat::Zero(3 * trigSize, 1);
            //    for (UInt i = 0; i != trigSize; ++i)
            //    {
            //        Vec3 grad = this->contactMeshBarrierGrad[i];
            //        analyticBarrierGrad.block<3, 1>(3 * i, 0) = grad;
            //    }
            //    DMat analyticBarrierGradT = analyticBarrierGrad.transpose();
            //    this->AppendToLog(ToStr(analyticBarrierGradT) + "\n");
            //}

            //{
            //    this->AppendToLog("===== Barrier Hessian Test [Trig] =====\n\n");

            //    for (UInt i = 0; i != numericBarrierHess_Trig.rows(); ++i)
            //        for (UInt j = 0; j != numericBarrierHess_Trig.cols(); ++j)
            //            if (std::abs(numericBarrierHess_Trig(i, j)) < 1e-2)
            //                numericBarrierHess_Trig(i, j) = 0;
            //    this->AppendToLog("Numerical Barrier Hess [Trig]: \n" + ToStr(numericBarrierHess_Trig) + "\n");

            //    SpMatData analyticBarData = this->contactMeshBarrierHess.toTriplets(this->contactMesh.boundaryTypes);
            //    SpMat analyticBarrierHess(3 * trigSize, 3 * trigSize);
            //    analyticBarrierHess.setFromTriplets(analyticBarData.begin(), analyticBarData.end());
            //    DMat analyticBarHessD = SparseToDense(analyticBarrierHess);
            //    this->AppendToLog("Analytic Barrier Hess [Trig]: \n" + ToStr(analyticBarHessD) + "\n");
            //}

            //DMat numericBarrierGrad_BS = this->NumericalBarrierGrad_BS();
            //DMat numericBarrierHess_BS = this->NumericalBarrierHess_BS();
            //{
            //    this->AppendToLog("===== Barrier Gradient Test [BS] =====\n\n");
            //    DMat numericBarGrad_BS_T = numericBarrierGrad_BS.transpose();
            //    this->AppendToLog("Numerical Barrier Grad [BS]: \n" + ToStr(numericBarGrad_BS_T) + "\n");

            //    UInt size = this->bsVertCntSum;
            //    DMat analyticBarrierGrad_BS = this->AccumulateTrigGrad(this->contactMeshBarrierGrad);
            //    DMat analyticBarrierGrad_BS_T = analyticBarrierGrad_BS.transpose();
            //    this->AppendToLog("Analytic Barrier Grad [BS]: \n" + ToStr(analyticBarrierGrad_BS_T) + "\n");
            //}

            //{
            //    this->AppendToLog("===== Barrier Hessian Test [BS] =====\n\n");
            //    for (UInt i = 0; i != numericBarrierHess_BS.rows(); ++i)
            //        for (UInt j = 0; j != numericBarrierHess_BS.cols(); ++j)
            //            if (std::abs(numericBarrierHess_BS(i, j)) < 1e-2)
            //                numericBarrierHess_BS(i, j) = 0;
            //    this->AppendToLog("Numerical Barrier Hess [BS]: \n" + ToStr(numericBarrierHess_BS) + "\n");

            //    const SpMatData barrierHessDataVec = this->contactMeshBarrierHess.toTriplets(this->contactMesh.boundaryTypes);
            //    tbb::concurrent_vector<SpMatEntry> barrierHessData;
            //    tbb::parallel_for(
            //        0, static_cast<Int>(barrierHessDataVec.size()), 1, [&](Int i)
            //        {
            //            const Eigen::Triplet<Float>& entry = barrierHessDataVec[i];
            //            barrierHessData.push_back(SpMatEntry(entry.row(), entry.col(), entry.value()));
            //        }
            //    );
            //    const UInt barrierHessEntryCnt = 81 * barrierHessData.size();
            //    SpMatData analyticBarrierHessData;
            //    analyticBarrierHessData.resize(barrierHessEntryCnt);
            //    this->FillBarrierHess(barrierHessData, analyticBarrierHessData, 0);

            //    UInt size = this->bsVertCntSum;
            //    SpMat analyticBarrierHess(3 * size, 3 * size);
            //    analyticBarrierHess.setFromTriplets(analyticBarrierHessData.begin(), analyticBarrierHessData.end());
            //    DMat analyticBarrierHessD = SparseToDense(analyticBarrierHess);

            //    this->AppendToLog("Analytic Barrier Hess [BS]: \n" + ToStr(analyticBarrierHessD) + "\n");
            //}
        }
    }

    DMat Solver::NumericalIPGrad(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos)
    {
        UInt size = bsVertCntSum;

        DMat grad = DMat::Zero(3 * size, 1);
        const Float h = 1e-8;

        std::vector<Vec3> perturbVertices = hypPos;

        for (UInt i = 0; i != size; ++i)
        {
            for (UInt index = 0; index != 3; ++index)
            {
                UInt nodeIndex = i;

                perturbVertices[i][index] += h;
                this->UpdateControlPoints(perturbVertices);
                Float upVal = this->IPVal(perturbVertices, estPos);

                perturbVertices[i][index] -= 2 * h;
                this->UpdateControlPoints(perturbVertices);
                Float downVal = this->IPVal(perturbVertices, estPos);

                perturbVertices[i][index] += h;
                this->UpdateControlPoints(perturbVertices);

                grad(3 * i + index, 0) = (upVal - downVal) / (2 * h);
            }
        }

        return grad;
    }

    DMat Solver::NumericalIPHess(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos)
    {
        UInt size = this->bsVertCntSum;

        DMat hess = DMat::Zero(3 * size, 3 * size);
        Float h = 1e-6;
        Float hSq = h * h;

        std::vector<Vec3> perturbVertices = hypPos;

        for (UInt i = 0; i != size; ++i)
        {
            for (UInt j = 0; j != size; ++j)
                for (UInt indexI = 0; indexI != 3; ++indexI)
                    for (UInt indexJ = 0; indexJ != 3; ++indexJ)
                    {
                        perturbVertices[i][indexI] += h;
                        perturbVertices[j][indexJ] += h;
                        this->UpdateControlPoints(perturbVertices);
                        Float upUpVal = this->IPVal(perturbVertices, estPos);

                        perturbVertices[j][indexJ] -= 2 * h;
                        this->UpdateControlPoints(perturbVertices);
                        Float upDownVal = this->IPVal(perturbVertices, estPos);

                        perturbVertices[i][indexI] -= 2 * h;
                        this->UpdateControlPoints(perturbVertices);
                        Float downDownVal = this->IPVal(perturbVertices, estPos);

                        perturbVertices[j][indexJ] += 2 * h;
                        this->UpdateControlPoints(perturbVertices);
                        Float downUpVal = this->IPVal(perturbVertices, estPos);


                        perturbVertices[i][indexI] += h;
                        perturbVertices[j][indexJ] -= h;
                        this->UpdateControlPoints(perturbVertices);

                        hess(3 * i + indexI, 3 * j + indexJ) = (upUpVal - upDownVal - downUpVal + downDownVal) / (4 * hSq);
                    }

            BSIPC_PEEK(i);
        }

        for (UInt i = 0; i != hess.cols(); ++i)
            for (UInt j = 0; j != hess.cols(); ++j)
                if (std::abs(hess(i, j)) < 1e-8)
                    hess(i, j) = 0;

        return hess;
    }
#endif

    void Solver::InitSimulation()
    {
#if defined BSIPC_USE_MKL
        BSIPC_WARN("Intel MKL using {} threads", mkl_get_max_threads());
#endif

        //BSIPC_INFO("Number of BSpline nodes: {}", this->bsVertCntSum);

        BSIPC_ASSERT(!this->targets.empty(), "No B-Spline Surface has been bound to the solver");

        // Initialize directory for writing stepCache
        BSIPC::OS::DeleteDir("steps");
        BSIPC::OS::CreateDir("steps");

        if (this->config.writeIPCMesh)
        {
            BSIPC::OS::DeleteDir("ipc_surf");
            BSIPC::OS::CreateDir("ipc_surf");
        }
        
        // if intended to use the quadrature point cached in `quad_cache`, change the directory name to `quad_data`
        if (this->config.cacheQuadPoints && !this->config.useCachedQuadPoints)
        {
            BSIPC::OS::DeleteDir("quad_cache");
            BSIPC::OS::CreateDir("quad_cache");
        }

        // Initialize the velocity in step cache
        this->stepCache.estVels.resize(this->bsVertCntSum, Vec3::Zero());

        // Initialize [BSTargetInfo]s (assume that [BSSurface]s has been bound to [BSTargetInfo]s)
        UInt bsVertIndexOffset = 0;
        UInt bsPatchIndexOffset = 0;
        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            BSTargetInfo& info = this->targets[idx];
            info.bsVertIndexOffset = bsVertIndexOffset;
            info.bsPatchIndexOffset = bsPatchIndexOffset;
            bsVertIndexOffset += info.target->GetControlPointCnt();
            bsPatchIndexOffset += info.target->GetPatchCnt();
        }

        this->localElasticityGrad.resize(this->bsPatchCntSum, Mat<27, 1>::Zero());
        this->localElasticityHess.resize(this->bsPatchCntSum, Mat<27, 27>::Zero());

        // Prepare the quadrature points, and material space related derivatives
        BSIPC_INFO("Start Initing Quadrature Points");
        for (UInt i = 0; i != this->targets.size(); ++i)
        {
            BSTargetInfo& info = this->targets[i];

            String quadDataPath = fmt::format("quad_data/surf_{}.yaml", info.index);
            String quadCachePath = fmt::format("quad_cache/surf_{}.yaml", info.index);

            if (this->config.useCachedQuadPoints)
                info.target->ReadQuadData(quadDataPath);
            else
            {
                BSIPC_PROFILE_SCOPE(
                    fmt::format("Quad Points for BS Surface [{}]", i),
                    info.target->PrepareQuadPoints();
                    );

                // This step needs to be done after quadrature points as we need integration 
                BSIPC_PROFILE_SCOPE(
                    "InitMassMatrix",
                    info.target->InitMassMatrix(config.density);
                );

                //info.target->TestQuadData(quadDataPath);
            }

            BSIPC_PEEK(info.target->GetMassMatEntries().size());
            BSIPC_PEEK(info.target->GetControlPointCnt());

            BSIPC_ASSERT(info.target->GetMassMatEntries().size() == info.target->GetControlPointCnt(), "mass matrix size sanity check");
            BSIPC_ASSERT(info.target->GetPatchStartIndices().size() == info.target->UDomainMax() * info.target->VDomainMax(), "patch size sanity check");

            if (this->config.cacheQuadPoints)
                info.target->WriteQuadCache(quadCachePath);
        }
        BSIPC_INFO("Finished Initing Quadrature Points");

        for (const BSTargetInfo& info : this->targets)
        {
            std::vector<Float> curMassMatEntries = info.target->GetMassMatEntries();
            this->massMatDiagEntries.insert(this->massMatDiagEntries.end(), curMassMatEntries.begin(), curMassMatEntries.end());
        }

        BSIPC_INFO("mass mat entry cnt: {}", this->massMatDiagEntries.size());

        if (this->config.continueOnStep)
            this->ReadFromStepLog();

        //this->AppendToLog("===== Diagonal Entries of Mass Matrix =====\n");
        //for (UInt i = 0; i != this->targets.size(); ++i)
        //{
        //    BSTargetInfo& info = this->targets[i];
        //    this->AppendToLog("BS Surface [" + std::to_string(i) + "]\n");

        //    for (Float entry : info.target->GetMassMatEntries())
        //        this->AppendToLog(std::to_string(entry) + "\n");
        //}
        //this->LogSep();

        for (auto& info : this->targets)
            info.target->UpdateQuadPointCache();

        this->InitMesh();
        this->InitHessCache();
        this->NormalizeSeamCoordinates();
        this->PrecalculateBendingHess();
        //this->FilterSeamContacts();
        this->PrecomputeKKTComponents();
        //this->PreAssembleElasticityHessTemplate();

        BSIPC_WARN("contact mesh vertex cnt: {}", this->contactMesh.vertices.size());

        if (this->config.writeIPCMesh)
            IPC::saveSurfaceMesh("ipc_surf/start.obj", this->contactMesh);

        for (const auto& info : this->targets)
        {
            BSIPC_ASSERT(info.Validate(), "Incomplete [BSTargetInfo]!");
        }

        this->WriteToStepLog("steps/start.yaml");
        BSIPC_INFO("Finished Initing Simulator");
    }

    void Solver::UpdateContactStiffness()
    {
        // enable hardcoding of contact stiffness
        if (this->config.contactStiffness.has_value())
            this->contactMesh.Kappa = this->config.contactStiffness.value();
        else
        {
            IPC::UpperBoundKappa(BSIPC_OUT this->contactStiffness, this->contactMesh.Hhat, this->contactMesh.bboxDiagSize2, this->contactMesh.meanMass);
            if (this->contactStiffness < 1e-16) {
                IPC::SuggestKappa(BSIPC_OUT this->contactStiffness, this->contactMesh.Hhat, this->contactMesh.bboxDiagSize2, this->contactMesh.meanMass);
            }
            IPC::InitKappa(this->contactMesh, this->ground, BSIPC_OUT this->contactStiffness);
            this->contactMesh.Kappa = this->contactStiffness;
        }

        BSIPC_PEEK(this->contactStiffness);
    }

    Bool Solver::StepForward()
    {
        BSIPC_PEEK(this->contactMesh.Hhat);

        Profiler::Start("Solve1");

        std::vector<Vec3> curTestBasisArray(this->bsVertCntSum);
        std::vector<Vec3> estPosArray(this->bsVertCntSum);
        std::vector<Vec3> prevPos(this->bsVertCntSum);

        for (UInt bsInd = 0; bsInd != this->targets.size(); ++bsInd)
        {
            const BSTargetInfo& info = this->targets[bsInd];
            const std::vector<Vec3>& curMeshPrevPos = info.target->GetControlPoints();

            tbb::parallel_for(
                0, static_cast<Int>(curMeshPrevPos.size()), 1, [&](Int i)
                {
                    UInt index = info.bsVertIndexOffset.value() + i;
                    prevPos[index] = curMeshPrevPos[i];
                    curTestBasisArray[index] = prevPos[index];
                }
            );
        }

        Float timestep = this->config.timestep;

        for (UInt bsInd = 0; bsInd != this->targets.size(); ++bsInd)
        {
            BSTargetInfo& info = this->targets[bsInd];

            UInt size = info.target->GetControlPointCnt();

            UInt startIndex = info.bsVertIndexOffset.value();

            tbb::parallel_for(
                0, static_cast<Int>(size), 1, [&](Int i)
                {
                    estPosArray[startIndex + i] = curTestBasisArray[startIndex + i] + timestep * this->stepCache.estVels[startIndex + i] + 
                        timestep * timestep * Vec3(0., 0., -this->config.gravConst);
                }
            );
        }

        if (this->animatedMeshCache.has_value() && this->stepCache.CurStep() < this->config.animatedMeshInfo.value().length)
        {
            const std::vector<Vec3>& prevFramePos = this->animatedMeshCache.value().curPos;
            const std::vector<Vec3>& curFramePos = this->animatedMeshPoses[this->stepCache.CurStep()].vertices;

            BSIPC_ASSERT(prevFramePos.size() == curFramePos.size(), "Inconsistent animated mesh vertex cnt");

            tbb::parallel_for(
                0, static_cast<Int>(this->animatedMeshCache.value().vertCnt), 1, [&](Int i)
                //for (UInt i = 0; i != this->animatedMeshCache.value().vertCnt; ++i)
                {
                    this->contactMesh.velocities[this->animatedMeshCache.value().startIdx + i] = (curFramePos[i] - prevFramePos[i]) / timestep;
                }
            );
        }

        BSIPC_PROFILE_SCOPE(
            "Core Solve",
            this->UpdateContactStiffness();

            // Semi-implicit lagged friction (Li et al. IPC 2020): snapshot tangent bases and lambdas
            // once per time step on the previous step's end-of-step geometry.  This is intentional —
            // do not move this call inside the Newton loop.  mesh.V_prev is not yet advanced (that
            // happens after the loop at line ~3012), so friction displacement correctly tracks
            // "displacement from end of step N-1 to current Newton iterate of step N".
            BSIPC_PROFILE_SCOPE(
                "Build Friction",
                IPC::BuildFriction(this->contactMesh, this->ground);
            )
            
            this->UpdateMovingDBC(curTestBasisArray, this->stepCache.CurStep());

            while (true)
            {
                IPC::PreSolveAction(this->contactMesh);

                //try
                //{
                    this->SolveSubIP(curTestBasisArray, estPosArray);
                //}
                //catch (const Exception& e)
                //{
                    //this->AppendToLog(std::string("Exception in [SolveSubIP]: \n") + e.Info() + "\n");
                    //BSIPC_ERROR("Exception in [SolveSubIP]: \n{}", e.Info());
                    ////return true;
                //}

                if (IPC::PostSolveAction(this->ground, this->contactMesh))
                    break;
            }

            this->UpdateControlPoints(curTestBasisArray);
        )

        Profiler::Check("Solve1");
        UInt curSimlTime = Profiler::GetTime("Solve1").value();
        Profiler::Stop("Solve1");

        this->LogSep();
        this->AppendToLog("Step: " + std::to_string(this->stepCache.CurStep()) +
            " Newton Iteration: " + std::to_string(this->stepCache.newtonIterCnt) + "\n");

        //Float maxZ = 0.;
        //const BSTargetInfo& info = this->targets[0];
        //for (Float u = 0; u < info.target->UDomainMax(); u += 0.1)
        //    for (Float v = 0; v < info.target->VDomainMax(); v += 0.1)
        //    {
        //        Vec3 pos = info.target->EvalAt(u, v);
        //        if (std::abs(pos[2]) > maxZ)
        //            maxZ = std::abs(pos[2]);
        //    }

        this->stepCache.NextStep(this->bsVertCntSum);

        /// Update Velocity & Position
        this->contactMesh.V_prev = this->contactMesh.vertices;
        IPC::ComputeXTilde(this->contactMesh);

        for (UInt bsInd = 0; bsInd != this->targets.size(); ++bsInd)
        {
            UInt size = this->targets[bsInd].target->GetControlPointCnt();

            tbb::parallel_for(
                0, static_cast<Int>(size), 1, [&](Int i)
                {
                    UInt index = this->targets[bsInd].bsVertIndexOffset.value() + i;
                    this->stepCache.estVels[index] = (curTestBasisArray[index] - prevPos[index]) / timestep;
                }
            );
        }

        /// Separate the steps
#if not defined BSIPC_PERF_TEST
        BSIPC_PROFILE_SCOPE(
            "Write to files",
            const String filename = "steps/step_" + std::to_string(this->stepCache.CurStep()) + ".yaml";
            this->WriteToStepLog(filename);
        )

        if (this->config.writeIPCMesh)
            IPC::saveSurfaceMesh("ipc_surf/surf_" + std::to_string(this->stepCache.CurStep()) + ".obj", this->contactMesh);
#endif

        BSIPC_PROFILE_SCOPE(
            "Propagate vel",
            this->PropagateVelToTrigMesh(stepCache.estVels);
        )

        /// Affliate updates
#if not defined BSIPC_DISABLE_RENDERER
        BSIPC_PROFILE_SCOPE(
            "UpdateRenderMesh",
            this->UpdateRenderMesh();
        )
#endif

        this->simlTime += curSimlTime;
        this->AppendToLog(fmt::format("## Step: {}, Current Step Time: {}, Summed Time: {}\n", 
            this->stepCache.CurStep(), curSimlTime, this->simlTime));

        BSIPC_WARN("Step: {}, Time: {}; Sum time: {}", this->stepCache.CurStep(), curSimlTime, this->simlTime);

        if (this->stepCache.CurStep() >= this->config.maxSteps)
            return true;
        else
            return false;
    }

    void Solver::SolveSubIP(BSIPC_OUT std::vector<Vec3>& curTestBasis, const std::vector<Vec3>& estPos)
    {
        Float timestep = this->config.timestep;
        DMat searchDir = DMat::Zero(3 * this->bsVertCntSum, 1);
        Float residue = 0.;

        Bool reachedToI = false;

        do
        {
            BSIPC_INFO_SEP

            Profiler::Start("single iter");

            BSIPC_PROFILE_SCOPE(
                "init control",
                // For vertices fixed by non-BS-surface mesh, we still need to provide velocity (set to 0)
                tbb::parallel_for(
                    0, static_cast<Int>(this->contactMeshSearchDir.size()), 1, [&](Int i)
                    {
                        this->contactMeshSearchDir[i] = Vec3::Zero();
                    }
                );
                this->UpdateControlPoints(curTestBasis);

               this->spatialHash.build(this->contactMesh, this->contactMesh.averageEdgeLength);
                IPC::BuildCollisionSets(this->contactMesh, this->spatialHash, this->ground, false);
            )

            this->includeAnimatedMeshInSystem = this->ActivateAnimatedMesh();
            //this->includeAnimatedMeshInSystem = true;
            std::vector<Vec3> curAnimatedMeshPrevPos;
            if (this->includeAnimatedMeshInSystem)
            {
                curAnimatedMeshPrevPos = this->animatedMeshCache.value().curPos;
                this->AppendToLog(">>>>> Including animation mesh in the system.");
            }

            BSIPC_WARN("step {}", this->stepCache.CurStep());

            BSIPC_PROFILE_SCOPE(
                "Barrier Grad/Hess",
                IPC::ComputeBarrierFrictionGradHess(this->contactMesh, this->contactMeshBarrierGrad, this->contactMeshBarrierHess, this->ground);
            )

            // Numeric test of barrier gradient/hessian
#if defined BSIPC_NUMERIC_TEST
            this->BarrierGradHessNumericTest();
#endif

            // P <- SPDProjection(Hess E)
            BSIPC_PROFILE_SCOPE(
                "IPHess",
                SpMat spdHess = IPHess(curTestBasis, estPos);
            )

            //BSIPC_TRACE(">> Finish computing Hessian");

            // p <- -P^{-1} * Grad E
            BSIPC_PROFILE_SCOPE(
                "IPGrad",
                DMat grad = IPGrad(curTestBasis, estPos);
            )

            ValidateSparseMatrix(spdHess);
            DMat searchDir;

            BSIPC_PROFILE_START("SolveLinearSystem");

            Bool success = true;
            if (this->seamConstraintCnt > 0)
                success = this->SolveKKTSystemWithDofReduction(spdHess, grad, BSIPC_OUT searchDir);
            else
            {
#if defined BSIPC_USE_EIGEN_SOLVER
                //BSIPC_INFO("hess size: {}x{}, grad size: {}", spdHess.rows(), spdHess.cols(), grad.rows());
                success = SolveLinearSystemLDLT(spdHess, -grad, searchDir);
                //DMat searchDir = SolveLinearSystemEigenCG(spdHess, -grad);
#elif defined BSIPC_USE_SUITESPARSE_SOLVER
                BSIPC_INFO("Solving using Suitesparse CHOLMOD");
                CholmodSolver cholmodSolver;
                DMat searchDir = cholmodSolver.SolveLinearSystem(spdHess, -grad);
                //DMat searchDir = SolveLinearSystemSuiteSparse(spdHess, -grad);
#else
#error "Need to specify which linear solver to use."
#endif
            }

            if (!success)
                this->AppendToLog("[ERROR] Linear Solver failed even with SimpliclalLDLT fallback!\n");

            BSIPC_PROFILE_END("SolveLinearSystem");

            BSIPC_WARN("grad.*searchDir: {}", DVec(grad).dot(DVec(searchDir)));
            this->AppendToLog(fmt::format("grad.*searchDir: {}\n", DVec(grad).dot(DVec(searchDir))));

            // Backtracking Filtering Line search
            BSIPC_PROFILE_SCOPE(
                "PropagateSearchDirToTrigMesh",
                this->PropagateSearchDirToTrigMesh(searchDir);
            )

            BSIPC_TRACE(">> Finish propagating search direction");

            // NOTE: EnvSelfCCDBound internally calls sh.build(mesh, moveDir, ...) — swept-CCD variant.
            // On return, spatialHash is in swept state. PreLSIntersectionCCDBound immediately rebuilds
            // the static hash, so this is self-correcting; do not reorder these two calls.
            BSIPC_PROFILE_SCOPE(
                "EnvCCD",
                Float alpha = IPC::EnvSelfCCDBound(this->contactMesh, this->spatialHash, this->ground, this->contactMeshSearchDir, this->contactStiffness);
            )

            std::vector<Vec3> curContactMeshVerts = this->contactMesh.vertices;
            Bool intersected;

            BSIPC_PROFILE_SCOPE(
               "Pre LS CCD",
               alpha = IPC::PreLSIntersectionCCDBound(
                   this->contactMesh, this->spatialHash, this->ground, alpha, curContactMeshVerts, this->contactMeshSearchDir, intersected
               );
            )
            if (intersected)
            {
                this->AppendToLog("## INTERSECTED!\n"); 
            }

            Float preLSAlpha = alpha;

            BSIPC_INFO("Start Search Length: {}", alpha);
            this->AppendToLog(fmt::format("Start Search Length: {}\n", alpha));

            //BSIPC_PROFILE_SCOPE(
                //"mid control",

                BSIPC_PROFILE_SCOPE(
                    "IPVal",
                    Float prevEnergy = IPVal(curTestBasis, estPos);
                )

                std::vector<Vec3> curTestPos(this->bsVertCntSum);
                tbb::parallel_for(
                    0, static_cast<Int>(this->bsVertCntSum), 1, [&](Int i)
                    {
                        curTestPos[i] = curTestBasis[i] + alpha * Vec3(searchDir(3 * i), searchDir(3 * i + 1), searchDir(3 * i + 2));
                    }
                );
                this->UpdateControlPoints(curTestPos);

                if (this->includeAnimatedMeshInSystem)
                {
                    //BSIPC_ASSERT(
                    //    searchDir.size() == 3 * (this->bsVertCntSum + this->animatedMeshCache.value().vertCnt + this->seamConstraintCnt), 
                    //    "[pre-Update] search direction size mismatch"
                    //);

                    tbb::parallel_for(
                        0, static_cast<Int>(this->animatedMeshCache.value().vertCnt), 1, [&](Int i)
                        {
                            UInt startIdx = 3 * (this->bsVertCntSum + i);
                            this->animatedMeshCache.value().curPos[i] = curAnimatedMeshPrevPos[i] + 
                                alpha * Vec3(searchDir(startIdx), searchDir(startIdx + 1), searchDir(startIdx + 2));
                        }
                    );
                    this->UpdateAnimatedMesh(this->animatedMeshCache.value().curPos);
                }
            //)

            this->stepCache.lineSearchIterCnt = 0;

            BSIPC_PROFILE_SCOPE(
                "Line Search",
                // Rebuild hash and active sets before each IPVal check so barrier energy is evaluated
                // against the correct contact pairs at the current trial position.
                // IPVal reads Self_ActiveSet/Environment_ActiveSet directly; without this rebuild it
                // would see the pairs from curTestBasis (start of iter), not curTestPos — causing the
                // line search to miss energy increases from newly-activated contact pairs.
                while (true)
                {
                    this->spatialHash.build(this->contactMesh, this->contactMesh.averageEdgeLength);
                    IPC::BuildCollisionSets(this->contactMesh, this->spatialHash, this->ground, false);

                    Float curEnergy = IPVal(curTestPos, estPos);
                    BSIPC_INFO("ls index: {}, cur: {:02.8f}, prev: {:02.8f}", this->stepCache.lineSearchIterCnt, curEnergy, prevEnergy);
                    this->AppendToLog(fmt::format("LS iter: {}, alpha: {:.8f}, cur: {:.8f}, prev: {:.8f}\n",
                        this->stepCache.lineSearchIterCnt, alpha, curEnergy, prevEnergy));

                    if (curEnergy <= prevEnergy)
                        break;

                    if (this->stepCache.lineSearchIterCnt >= this->config.maxLSIterations)
                    {
                        BSIPC_WARN("Line search: reached max iterations ({})", this->config.maxLSIterations);
                        this->AppendToLog(">>>> Line search: reached max iterations\n");
                        break;
                    }

                    ++this->stepCache.lineSearchIterCnt;
                    alpha *= 0.5;

                    tbb::parallel_for(
                        0, static_cast<Int>(this->bsVertCntSum), 1, [&](Int i)
                        {
                            curTestPos[i] = curTestBasis[i] + alpha * Vec3(searchDir(3 * i), searchDir(3 * i + 1), searchDir(3 * i + 2));
                        }
                    );
                    this->UpdateControlPoints(curTestPos);

                    if (this->includeAnimatedMeshInSystem)
                    {
                        tbb::parallel_for(
                            0, static_cast<Int>(this->animatedMeshCache.value().vertCnt), 1, [&](Int i)
                            {
                                UInt startIdx = 3 * (this->bsVertCntSum + i);
                                this->animatedMeshCache.value().curPos[i] = curAnimatedMeshPrevPos[i] +
                                    alpha * Vec3(searchDir(startIdx), searchDir(startIdx + 1), searchDir(startIdx + 2));
                            }
                        );
                        this->UpdateAnimatedMesh(this->animatedMeshCache.value().curPos);
                    }
                }
            )

            BSIPC_PROFILE_SCOPE(
                "Post LS CCD",
                Float postLSAlpha = IPC::PostLSIntersectionCCDBound(
                   preLSAlpha, alpha, this->contactMeshSearchDir, this->contactMesh, this->spatialHash, this->ground, curContactMeshVerts
                );

                if (postLSAlpha != alpha)
                {
                   alpha = postLSAlpha;

                   curTestPos.clear();
                   curTestPos.resize(bsVertCntSum);
                   
                   tbb::parallel_for(
                       0, static_cast<Int>(this->bsVertCntSum), 1, [&](Int ind)
                       {
                           curTestPos[ind] = curTestBasis[ind] + alpha * Vec3(searchDir(3 * ind), searchDir(3 * ind + 1), searchDir(3 * ind + 2));
                       }
                   );

                   if (this->includeAnimatedMeshInSystem)
                   {
                       BSIPC_ASSERT(searchDir.size() == 3 * (this->bsVertCntSum + this->animatedMeshCache.value().vertCnt), "search direction size mismatch");

                       tbb::parallel_for(
                           0, static_cast<Int>(this->animatedMeshCache.value().vertCnt), 1, [&](Int i)
                           {
                               UInt startIdx = 3 * (this->bsVertCntSum + i);
                               this->animatedMeshCache.value().curPos[i] = curAnimatedMeshPrevPos[i] +
                                   alpha * Vec3(searchDir(startIdx), searchDir(startIdx + 1), searchDir(startIdx + 2));
                           }
                       );
                       this->UpdateAnimatedMesh(this->animatedMeshCache.value().curPos);
                   }
                }

                curTestBasis = curTestPos;

                // PostLSIntersectionCCDBound calls IPC::stepForward when alpha < preLSAlpha, which
                // advances contactMesh.vertices via a linear trig-mesh step rather than spline
                // evaluation.  postLSAlpha always equals alpha (the CCD-reduction loop inside is
                // commented out), so testing postLSAlpha != alpha would never fire.  Instead test
                // alpha < preLSAlpha — the condition under which stepForward actually ran — and
                // re-sync the spline state before PostLineSearch snapshots closeConstraintID/Val.
                if (alpha < preLSAlpha)
                {
                    this->UpdateControlPoints(curTestPos);
                    this->spatialHash.build(this->contactMesh, this->contactMesh.averageEdgeLength);
                    IPC::BuildCollisionSets(this->contactMesh, this->spatialHash, this->ground, false);
                }
            )

            BSIPC_PROFILE_SCOPE(
                "PostLineSearch",
                IPC::PostLineSearch(this->contactMesh, this->ground, alpha, this->contactStiffness);
            )

            BSIPC_PROFILE_SCOPE(
                "post misc",

                ++this->stepCache.newtonIterCnt;
                UInt sysDoF = 3 * this->bsVertCntSum;
                if (this->includeAnimatedMeshInSystem)
                    sysDoF += 3 * this->animatedMeshCache.value().vertCnt;
                residue = searchDir.topRows(sysDoF).lpNorm<Eigen::Infinity>() / timestep;

                BSIPC_INFO("residue: {}", residue);

                this->AppendToLog("Line Search Iteration: " + std::to_string(this->stepCache.lineSearchIterCnt) +
                    " step scale (alpha): " + std::to_string(alpha) +
                    " residue: " + std::to_string(residue) +
                    "\n");
                this->AppendToLog("-------------------------------\n");
            )

            if (this->config.writeIterMesh)
            {
                String filename = fmt::format("steps/step_{}_iter_{}.yaml", this->stepCache.CurStep(), this->stepCache.newtonIterCnt);
                this->WriteToStepLog(filename);
                if (this->config.writeIPCMesh)
                    IPC::saveSurfaceMesh(fmt::format("ipc_surf/surf_{}_iter_{}.obj", this->stepCache.CurStep(), this->stepCache.newtonIterCnt), this->contactMesh);
            }

            Profiler::Check("single iter");
            BSIPC_TRACE("{}: {}", "single iter", Profiler::GetTime("single iter").value());
            Profiler::Stop("single iter");

            if (this->stepCache.newtonIterCnt >= this->config.maxNewtonIterations)
            {
                BSIPC_WARN("Reach max Newton iterations: {}", this->config.maxNewtonIterations);
                this->AppendToLog(">>>> Reach max Newton iterations: " + std::to_string(this->config.maxNewtonIterations) + "\n");
                break;
            }

        } while (residue > this->config.tolerance);
    }

    void Solver::SolveKKTSystemSchur(const SpMat& hess, const DMat& rhs, BSIPC_OUT DMat& x) const
    {
#if defined BSIPC_USE_MKL
        const SpMat& sparseC = this->includeAnimatedMeshInSystem ? this->seamConstraintCoeffsExtended : this->seamConstraintCoeffs;
        BSIPC_ASSERT(hess.rows() == hess.cols(), "Hessian should be square");
        UInt hessSize = hess.rows();
        UInt dofSize = hessSize - 3 * this->seamConstraintCnt;

        SpMatData selectionData;
        selectionData.reserve(dofSize);
        for (UInt i = 0; i != dofSize; ++i)
            selectionData.push_back(SpMatEntry(i, i, 1.));
        SpMat selection(dofSize, hessSize);
        selection.setFromTriplets(selectionData.begin(), selectionData.end());

        SpMat H(dofSize, dofSize);
        H = selection * hess * selection.transpose();
        DMat g = selection * rhs;

        Eigen::PardisoLDLT<SpMat> ldlt_H;
        ldlt_H.compute(H);
        BSIPC_ASSERT(ldlt_H.info() == Eigen::Success, "ldlt_H.compute failed");

        DMat C_T = sparseC.transpose();
        DMat S = -sparseC * ldlt_H.solve(C_T);

        Eigen::PartialPivLU<DMat> lu_S;
        lu_S.compute(S);

        DMat H_inv_grad = ldlt_H.solve(g);
        DMat lambda = lu_S.solve(-sparseC * H_inv_grad);
        DMat xRhs = g - sparseC.transpose() * lambda;

        x = ldlt_H.solve(xRhs);

        DMat diff = H * x - g;
        BSIPC_PEEK(diff.lpNorm<Eigen::Infinity>());

        if (diff.lpNorm<Eigen::Infinity>() > 0.5 * g.lpNorm<Eigen::Infinity>())
        {
            BSIPC_CRITICAL("Schur PardisoLDLT gives wrong result!");
        }
#else
        BSIPC_WARN("MKL disabled. Falling back to SimplicialLDLT solving");
        Bool success = SolveLinearSystemLDLT(hess, rhs, x);
#endif
    }

    Bool Solver::SolveKKTSystemWithDofReduction(const SpMat& hess, const DMat& grad, BSIPC_OUT DMat& x) const
    {
        const SpMat& R = this->includeAnimatedMeshInSystem ? this->rowReorderingMatExtended : this->rowReorderingMat;
        const SpMat& extendedZ = this->includeAnimatedMeshInSystem ? this->constrainedDofProjMatExtended : this->constrainedDofProjMat;
        BSIPC_INFO("R: {}x{}, extendedZ: {}x{}, hess: {}x{}", R.rows(), R.cols(), extendedZ.rows(), extendedZ.cols(), hess.rows(), hess.cols());

        SpMat projectedH;
        BSIPC_PROFILE_START("project Hessian MKL");
        ProjectKKTHessian(extendedZ, R, hess, projectedH);
        BSIPC_PROFILE_END("project Hessian MKL");

        BSIPC_PROFILE_START("project gradient");
        const DMat projectedGradNegated = -extendedZ.transpose() * R * grad;
        BSIPC_PROFILE_END("project gradient");
        DMat reducedSearchDir;
        Bool success = SolveLinearSystemLDLT(projectedH, projectedGradNegated, BSIPC_OUT reducedSearchDir);
        BSIPC_PROFILE_START("project searchDir");
        x = R.transpose() * extendedZ * reducedSearchDir;
        BSIPC_PROFILE_END("project searchDir");

        return success;
    }

    void Solver::PropagateSearchDirToTrigMesh(const DMat& searchDir)
    {
        UInt searchDirSize = this->bsVertCntSum;
        if (this->includeAnimatedMeshInSystem)
            searchDirSize += this->animatedMeshCache.value().vertCnt;

        //BSIPC_ASSERT(searchDir.rows() == 3 * searchDirSize + 3 * this->seamConstraintCnt, "Gradient size mismatch");

        tbb::parallel_for(
            0, static_cast<Int>(this->trigVertCntSum), 1, [&](Int i)
            {
                Vec3 localGrad = Vec3::Zero();

                // Here need to invert the direction as in IPC search direction is defined in the opposite direction
                for (const SourceInfo& source : this->linearVertexSources[i])
                    localGrad -= Vec3(searchDir.block<3, 1>(3 * source.index, 0)) * source.coeff;

                this->contactMeshSearchDir[i] = localGrad;
            }
        );

        if (this->includeAnimatedMeshInSystem)
        {
            tbb::parallel_for(
                0, static_cast<Int>(this->animatedMeshCache.value().vertCnt), 1, [&](Int i)
            //for (UInt i = 0; i != this->animatedMeshCache.value().vertCnt; ++i)
                {
                    this->contactMeshSearchDir[this->animatedMeshCache.value().startIdx + i] = -searchDir.block<3, 1>(3 * (this->bsVertCntSum + i), 0);
                }
            );
        }
        else
        {
            if (this->animatedMeshCache.has_value())
            {
                // Set the search direction for animated mesh vertices to zero
                tbb::parallel_for(
                    0, static_cast<Int>(this->animatedMeshCache.value().vertCnt), 1, [&](Int i)
                //for (UInt i = 0; i != this->animatedMeshCache.value().vertCnt; ++i)
                    {
                        this->contactMeshSearchDir[this->animatedMeshCache.value().startIdx + i] = Vec3::Zero();
                    }
                );
            }
        }
    }

    void Solver::PropagateVelToTrigMesh(const std::vector<Vec3>& vel)
    {
        UInt trigVerticesCnt = this->trigVertCntSum;

        tbb::parallel_for(
            0, static_cast<Int>(trigVerticesCnt), 1, [&](Int i)
            {
                this->contactMesh.velocities[i] = Vec3::Zero();

                for (UInt j = 0; j != this->linearVertexSources[i].size(); ++j)
                {
                    const SourceInfo& source = this->linearVertexSources[i][j];
                    this->contactMesh.velocities[i] += source.coeff * vel[source.index];
                }
            }
        );
    }

    DMat Solver::AccumulateTrigGrad(const std::vector<Vec3>& trigGrad) const
    {
        // TODO test which one is faster: 
        //  - Fill in parallel with [SpMatData] and then generate a dense view
        //  - Use non-parallel for and fill ordinarily

        BSIPC_ASSERT(
            trigGrad.size() == this->contactMesh.vertices.size(), "Barrier gradient vector should have the same size as mesh vertices"
        );

        SpMatData barrierGradData;
        UInt size = this->bsVertCntSum;
        UInt verticesCnt = this->trigVertCntSum;
        const UInt sourceCnt = 9;
        barrierGradData.resize(sourceCnt * verticesCnt * 3); // default constructor init this to (0, 0, 0)

        tbb::parallel_for(
            0, static_cast<Int>(verticesCnt), 1, [&](Int i)
            {
                Vec3 localGrad = trigGrad[i];
                for (UInt sourceIndex = 0; sourceIndex != this->linearVertexSources[i].size(); ++sourceIndex)  // this->vertexSources[i].size() = 9
                {
                    const SourceInfo& source = this->linearVertexSources[i][sourceIndex];
                    Vec3 sourceGrad = source.coeff * localGrad;
                    for (UInt j = 0; j != 3; ++j)
                        barrierGradData[27 * i + 3 * sourceIndex + j] = SpMatEntry(3 * source.index + j, 0, sourceGrad[j]);
                }
            }
        );

        SpMat barrierGradSparse(3 * size, 1);
        barrierGradSparse.setFromTriplets(barrierGradData.begin(), barrierGradData.end());

        // Invert the sign here
        DMat barrierGradDense = SparseToDense(barrierGradSparse);
        return barrierGradDense;
    }
    
    void Solver::AppendToLog(String str) const
    {
        std::ofstream ofs;
        ofs.open("log.txt", std::ios::app | std::ios::out);
        ofs << str;
        ofs.close();
    }

    ///  Helpers
    void Solver::WriteToLog() const
    {
        std::ofstream ofs;
        ofs.open("log.txt", std::ios::app | std::ios::out);

        //ofs << "======================== Step ========================\n";

        //ofs << "Control Points:\n";
        //auto& vertices = this->target->GetControlPoints();
        //for (auto& vertex : vertices)
        //    ofs << vertex[0] << " " << vertex[1] << " " << vertex[2] << " " << "\n";

        //ofs << "== Mesh:\n" << "==== Vertices:\n";
        //auto meshVertices = this->mesh.GetVertices();
        //auto meshIndices = this->mesh.GetIndices();

        //for (auto& vertex : meshVertices)
        //    ofs << vertex.pos.x << " " << vertex.pos.y << " " << vertex.pos.z << "\n";

        //ofs << "==== Indices:\n";
        //for (auto& index : meshIndices)
        //    ofs << index << "\n";

        ofs.close();
    }

    void Solver::WriteToStepLog(const String filename) const
    {
        fkyaml::node root = this->config.configRootNode;

        // Log BSSurface control points (in the world space)
        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            fkyaml::node& curSurfNode = root[Config::BS_SURFACES_STR][idx];
            std::vector<fkyaml::node> curControlNodes;

            const std::vector<Vec3>& controlPoints = this->targets[idx].target->GetControlPoints();
            for (auto iter = controlPoints.cbegin(); iter != controlPoints.cend(); ++iter)
                curControlNodes.push_back(fkyaml::node{ (*iter)[0], (*iter)[1], (*iter)[2] });
            curSurfNode[Config::CONTROL_NODES_STR] = curControlNodes;

            UInt bsOffset = this->targets[idx].bsVertIndexOffset.value();
            UInt bsSize = this->targets[idx].target->GetControlPointCnt();
            std::vector<fkyaml::node> curEstVelNodes;
            const std::vector<Vec3>& estVels = this->stepCache.estVels;
            for (UInt i = bsOffset; i != bsOffset + bsSize; ++i)
                curEstVelNodes.push_back(fkyaml::node{ estVels[i][0], estVels[i][1], estVels[i][2] });
            curSurfNode[Config::EST_VELOCITY_STR] = curEstVelNodes;
        }

        // Log animated meshes, if there is any
        if (this->animatedMeshCache.has_value())
        {
            fkyaml::node& animatedMeshNode = root[Config::ANIMATED_MESHES_STR];
            std::vector<fkyaml::node> curAnimatedMeshNodes;
            const std::vector<Vec3>& animatedMeshPos = this->animatedMeshCache.value().curPos;

            for (auto iter = animatedMeshPos.cbegin(); iter != animatedMeshPos.cend(); ++iter)
                curAnimatedMeshNodes.push_back(fkyaml::node{ (*iter)[0], (*iter)[1], (*iter)[2] });

            animatedMeshNode[Config::VERTICES_STR] = curAnimatedMeshNodes;
        }

        // Log StepCache: current step count, and estimated velocity
        UInt curStep = this->stepCache.CurStep();
        root[Config::CURRENT_STEP_STR] = curStep;

        // Write to file
        std::ofstream ofs(filename);
        if (ofs)
            ofs << root;
    }

    void Solver::ReadFromStepLog(const String path)
    {
        std::ifstream ifs(path);

        const SolverConfig& config = this->config;
        BSIPC_ASSERT(
            this->targets.size() == config.surfacesInfo.size(),
            "Solver should have the same number of surfaces as specified in config after BSSurfaceMgr parsing"
        );

        try
        {
            if (!ifs)
            {
                std::string msg = "error opening config file " + path;
                throw fkyaml::exception(msg.c_str());
            }
            fkyaml::node root = fkyaml::node::deserialize(ifs);

            // Read StepCache
            UInt startStep = 0;
            LoadDataFromYaml(startStep, root, Config::CURRENT_STEP_STR);
            this->stepCache.StartFromStep(startStep);
            this->stepCache.estVels.clear();

            // Read control points and velocity
            std::vector<Vec3> globalControlPts;
            for (UInt idx = 0; idx != this->targets.size(); ++idx)
            {
                const fkyaml::node& curSurfNode = root[Config::BS_SURFACES_STR][idx];
                const fkyaml::node& controlPtNode = curSurfNode[Config::CONTROL_NODES_STR];
                const fkyaml::node& estVelNode = curSurfNode[Config::EST_VELOCITY_STR];

                BSIPC_ASSERT(
                    controlPtNode.is_sequence() || controlPtNode.is_string(), 
                    "\"control_nodes\" node should either be a sequence (of Vec3) or a string (filename, with target file containing a sequence of Vec3)"
                );

                for (const fkyaml::node& node : controlPtNode)
                {
                    BSIPC_ASSERT(node.is_sequence() && node.size() == 3, "Control point node is not a sequence or has size different from 3");
                    globalControlPts.emplace_back(node[0].get_value<Float>(), node[1].get_value<Float>(), node[2].get_value<Float>());
                }

                for (const fkyaml::node& node : estVelNode)
                {
                    BSIPC_ASSERT(node.is_sequence() && node.size() == 3, "Velocity node is not a sequence or has size different from 3");
                    this->stepCache.estVels.emplace_back(node[0].get_value<Float>(), node[1].get_value<Float>(), node[2].get_value<Float>());
                }
            }

            BSIPC_ASSERT(
                globalControlPts.size() == this->bsVertCntSum,
                fmt::format("Control points size mismatch: config({}), step log({})", globalControlPts.size(), this->bsVertCntSum)
            );

            BSIPC_ASSERT(this->stepCache.estVels.size() == this->bsVertCntSum, "Velocity size mismatch");

            this->UpdateControlPoints(globalControlPts);

            // Read animated mesh vertices, if any
            if (root.contains(Config::ANIMATED_MESHES_STR))
            {
                BSIPC_WARN("Reading animated mesh info from config!");

                const fkyaml::node& animatedMeshNode = root[Config::ANIMATED_MESHES_STR];
                const fkyaml::node& animatedMeshVertNode = animatedMeshNode[Config::VERTICES_STR];
                BSIPC_ASSERT(
                    animatedMeshVertNode.is_sequence() || animatedMeshVertNode.is_string(),
                    "\"vertices\" node should either be a sequence (of Vec3) or a string (filename, with target file containing a sequence of Vec3)"
                );
                std::vector<Vec3> animatedMeshPos;
                for (const fkyaml::node& node : animatedMeshVertNode)
                {
                    BSIPC_ASSERT(node.is_sequence() && node.size() == 3, "Animated mesh vertex node is not a sequence or has size different from 3");
                    animatedMeshPos.emplace_back(node[0].get_value<Float>(), node[1].get_value<Float>(), node[2].get_value<Float>());
                }

                this->animatedMeshCache = AnimatedMeshCache{
                    .startIdx = std::numeric_limits<UInt>::max(),
                    .vertCnt = std::numeric_limits<UInt>::max(),
                    .curPos = std::move(animatedMeshPos)
                };
            }
        }
        catch (fkyaml::exception& e)
        {
            BSIPC_WARN(e.what());
            BSIPC_WARN("Error loading step infos upon continuing simulation");
        }
    }

    void Solver::ReadFromStepLog()
    {
        const String filename = "config.yaml";
        this->ReadFromStepLog(filename);
    }

    void Solver::InitMesh()
    {
        std::vector<std::pair<UInt, UInt>> contactBSTrigSizes(this->targets.size());
        std::vector<std::pair<UInt, UInt>> renderBSTrigSizes(this->targets.size());

        // For contact on triangular mesh, the resolution of tessellation for contact is the same as the resolution of the mesh for rendering
        for (UInt i = 0; i != this->targets.size(); ++i)
        {
            const BSTargetInfo& info = this->targets[i];
            UInt uSize = static_cast<UInt>(info.target->UDomainMax() * config.perPatchMesh + 1);
            UInt vSize = static_cast<UInt>(info.target->VDomainMax() * config.perPatchMesh + 1);

            renderBSTrigSizes[i] = std::make_pair(uSize, vSize);
            contactBSTrigSizes[i] = std::make_pair(info.target->GetLinearResolutionX(), info.target->GetLinearResolutionY());
        }

        /** Indices of the vertices
         *  (vSize - 1) * uSize + 0  (vSize - 1) * uSize + 1  (vSize - 1) * uSize + 2  ...  (vSize - 1) * uSize + uSize - 1
         *  .                        .                        .                        ...  .
         *  .                        .                        .                        ...  .
         *  .                        .                        .                        ...  .
         *  0                        1                        2                        ...  uSize - 1
         */

        this->InitContactMesh(contactBSTrigSizes);
        
        // Reserve slots for vertices in fixed mesh
        this->contactMeshSearchDir.resize(this->contactMesh.vertices.size(), Vec3::Zero());

#if not defined BSIPC_DISABLE_RENDERER
        this->InitRenderMesh(renderBSTrigSizes);
        this->renderer->UpdateMesh(renderMesh);

        BSIPC_INFO("render mesh indices array size: {}", this->renderMesh.GetIndices().size());
        BSIPC_INFO("render mesh vertices cnt: {}", this->renderMesh.GetVertices().size());
#endif
    }

    void Solver::InitAnimatedMesh()
    {
        const AnimatedMeshInfo& info = this->config.animatedMeshInfo.value();

        // Pre-load all poses of the animated mesh
        this->animatedMeshPoses.resize(info.length);
        size_t pos = info.pattern.find("@");
        BSIPC_ASSERT(pos != String::npos, "Pattern for animated mesh must contain '@' to be replaced by frame index");
        for (UInt i = 0; i != info.length; ++i)
        {
            String filename = info.pattern;
            filename.replace(pos, 1, std::to_string(i));
            String path = fmt::format("models/{}/{}.obj", info.dirName, filename);

            this->animatedMeshPoses[i].load_triangleMesh(path, 1.0, Vec3(0.0, 0.0, 0.0));

            for (const LinearTransformation& transformation : info.transformations)
                for (Vec3& vertex : this->animatedMeshPoses[i].vertices)
                    transformation.ApplyTo(vertex);
        }

        BSIPC_WARN("CurStep: {}", this->stepCache.CurStep());

        // Stack the animated mesh vertices onto the contact mesh
        UInt curMeshVertCnt = this->contactMesh.vertices.size();
        UInt startAnimatedMeshIdx = 0;
        if (this->stepCache.CurStep() != 0)
            startAnimatedMeshIdx = std::min(this->stepCache.CurStep(), static_cast<UInt>(this->animatedMeshPoses.size())) - 1;
        const IPC::mesh3D& startMesh = this->animatedMeshPoses[startAnimatedMeshIdx];
        //const IPC::mesh3D& startMesh = this->animatedMeshPoses[0];
        for (UInt i = 0; i != startMesh.vertices.size(); ++i)
        {
            this->contactMesh.vertices.emplace_back(startMesh.vertices[i]);
            this->contactMesh.velocities.emplace_back(Vec3::Zero());
            this->contactMesh.Constraints.emplace_back(Mat<3, 3>::Identity());
            this->contactMesh.masses.emplace_back(0.);
            this->contactMesh.boundaryTypes.emplace_back(0);
        }

        for (const IVec3& trig : this->animatedMeshPoses[0].triangles)
        {
            IVec3 newTrig = trig + IVec3(curMeshVertCnt, curMeshVertCnt, curMeshVertCnt);
            this->contactMesh.triangles.emplace_back(newTrig);
        }

        UInt endMeshVertCnt = this->contactMesh.vertices.size();
        UInt animatedMeshStartIdx = curMeshVertCnt;
        UInt animatedMeshVertCnt = endMeshVertCnt - curMeshVertCnt;

        if (this->animatedMeshCache.has_value())    // if the position of the animated mesh has been specified in the step yaml log
        {
            AnimatedMeshCache& cache = this->animatedMeshCache.value();
            cache.startIdx = animatedMeshStartIdx;
            cache.vertCnt = animatedMeshVertCnt;

            BSIPC_ASSERT(
                cache.curPos.size() == cache.vertCnt,
                "Animated mesh cache size mismatch upon initialization"
            );

            this->UpdateAnimatedMesh(cache.curPos);
        }
        else
        {
            AnimatedMeshCache cache;
            cache.startIdx = animatedMeshStartIdx;
            cache.vertCnt = animatedMeshVertCnt;
            this->animatedMeshCache = cache;
            this->animatedMeshCache.value().curPos.resize(cache.vertCnt);

            tbb::parallel_for(
                0, static_cast<Int>(cache.vertCnt), 1, [&](Int i)
                {
                    this->animatedMeshCache.value().curPos[i] = this->contactMesh.vertices[cache.startIdx + i];
                }
            );
        }
    }

    void Solver::InitContactMesh(std::vector<std::pair<UInt, UInt>> bsTrigSizes)
    {
        // Layout
        //  BS surface | animated mesh | fixed meshes

        // generate linear triangular mesh for B-spline surfaces
        this->bsVertexSupports.resize(this->bsVertCntSum);
        UInt patchCnt = this->targets.size();

        IPC::mesh3D& mesh = this->contactMesh;
        mesh.characterAnimation = this->config.disableInterpatchContact;
        //mesh.characterAnimation = true;
        mesh.vertices.clear();
        mesh.triangles.clear();
        mesh.patchVertexSep.resize(this->targets.size() + 1);
        mesh.patchEdgeSep.resize(this->targets.size() + 1);
        mesh.patchSurfaceSep.resize(this->targets.size() + 1);

        UInt contactTrigVertIndexOffset = 0;
        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            mesh.patchVertexSep[idx] = mesh.vertices.size();
            mesh.patchSurfaceSep[idx] = mesh.triangles.size();

            UInt trigVertCnt = 0;

            BSTargetInfo& info = this->targets[idx];    
            info.contactTrigVertIndexOffset = contactTrigVertIndexOffset;
            UInt uSize = bsTrigSizes[idx].first, vSize = bsTrigSizes[idx].second;
            UInt bsOffset = info.bsVertIndexOffset.value();

            BSIPC_INFO("usize: {}, vsize: {}", uSize, vSize);

            for (UInt j = 0; j != vSize; ++j)
            {
                for (UInt i = 0; i != uSize; ++i)
                {
                    Float u = static_cast<Float>(i) / static_cast<Float>(uSize - 1) * info.target->UDomainMax();
                    Float v = static_cast<Float>(j) / static_cast<Float>(vSize - 1) * info.target->VDomainMax();

                    Vec3 pos = info.target->EvalAt(u, v);

                    std::vector<SourceInfo> sources;

                    IVec2 supportStart = info.target->SupportStartOf(u, v);

                    Float coeff = 0;

                    for (UInt ctrlI = supportStart[0]; ctrlI < supportStart[0] + 3U && ctrlI < info.target->GetResolutionX(); ++ctrlI)
                    {
                        for (UInt ctrlJ = supportStart[1]; ctrlJ < supportStart[1] + 3U && ctrlJ < info.target->GetResolutionY(); ++ctrlJ)
                        {
                            UInt globalIndex = info.target->GetControlPointIndex(ctrlI, ctrlJ) + bsOffset;
                            Float coeffAtIJ = info.target->EvalCoeffAt(ctrlI, ctrlJ, u, v);
                            coeff += coeffAtIJ;
                            this->bsVertexSupports[globalIndex].emplace_back(coeffAtIJ, mesh.vertices.size());
                            sources.emplace_back(coeffAtIJ, globalIndex);
                        }
                    }

                    // ensures gradient/hessian filling process later is coherent with reserved size
                    while (sources.size() < 9)
                        sources.emplace_back(0, 0);

                    BSIPC_ASSERT(
                        sources.size() == 9, "Each point on B-spline surface must have 9 sources"
                    );

                    // ensure that the indices of vertices in mesh and sources align
                    this->linearVertexSources.emplace_back(sources);

                    ++trigVertCnt;
                    ++contactTrigVertIndexOffset;

                    mesh.vertices.emplace_back(pos);
                    mesh.velocities.emplace_back(Vec3::Zero());
                    mesh.Constraints.emplace_back(Mat<3, 3>::Identity());
                    mesh.masses.emplace_back(0.);
                    mesh.boundaryTypes.emplace_back(0);
                }
            }

            info.trigVertCnt = trigVertCnt;

            UInt trigOffset = info.contactTrigVertIndexOffset.value();
            for (UInt j = 0; j != vSize - 1; ++j)
            {
                for (UInt i = 0; i != uSize - 1; ++i)
                {
                    UInt lowerLeft = trigOffset + j * uSize + i;
                    UInt lowerRight = lowerLeft + 1;
                    UInt upperLeft = lowerLeft + uSize;
                    UInt upperRight = upperLeft + 1;

                    if ((i + j) % 2 == 0)
                    {
                        this->FormTrigInIPCMesh(lowerLeft, lowerRight, upperLeft);
                        this->FormTrigInIPCMesh(lowerRight, upperRight, upperLeft);
                    }
                    else
                    {
                        this->FormTrigInIPCMesh(lowerLeft, lowerRight, upperRight);
                        this->FormTrigInIPCMesh(lowerLeft, upperRight, upperLeft);
                    }
                }
            }
        }

        mesh.patchVertexSep[patchCnt] = mesh.vertices.size();
        mesh.patchSurfaceSep[patchCnt] = mesh.triangles.size();

        this->trigVertCntSum = static_cast<UInt>(mesh.vertices.size());

        this->contactMesh.friction = this->config.frictionMu;

        if (this->config.animatedMeshInfo.has_value())
            this->InitAnimatedMesh();

        // load all fixed linear meshes
        for (const FixedMeshInfo& info : this->config.fixedMeshInfos)
        {
            UInt startIndex = static_cast<UInt>(this->contactMesh.vertices.size());
            if (info.type == FixedMeshType::TETRAHEDRAL)
            {
                String meshPath = fmt::format("models/{}.msh", info.meshName);
                this->contactMesh.load_tetrahedraMesh(meshPath, 1.0, Vec3(0, 0, 0.0));
            }
            else if (info.type == FixedMeshType::TRIANGLE)
            {
                String meshPath = fmt::format("models/{}.obj", info.meshName);
                this->contactMesh.load_triangleMesh(meshPath, 1.0, Vec3(0, 0, 0.0));
            }
            //this->contactMesh.load_tetrahedraMesh("models/bunny.msh", 0.3, Vec3(0, 0, 0.0));
            UInt endIndex = static_cast<UInt>(this->contactMesh.vertices.size());

            UInt transformationCnt = info.transformations.size();

            for (const LinearTransformation& transformation : info.transformations)
                for (UInt ind = startIndex; ind != endIndex; ++ind)
                    transformation.ApplyTo(this->contactMesh.vertices[ind]);

            for (UInt ind = startIndex; ind != endIndex; ++ind)
            {
                this->contactMesh.boundaryTypes[ind] = 1;
                // if (info.type == FixedMeshType::TRIANGLE)
                // {
                    // this->contactMesh.boundaryTypes[ind] = 2;
                this->contactMesh.Constraints[ind].setZero();
                // }
            }
        }

        IPC::InitSettings(this->contactMesh);
        this->contactMesh.Hhat = this->config.contactThreshold * this->config.contactThreshold;

        BSIPC_WARN("active set size: {}", this->contactMesh.Self_ActiveSet.size());
        IPC::BuildModel(mesh, spatialHash, ground);
        BSIPC_WARN("active set size: {}", this->contactMesh.Self_ActiveSet.size());

        BSIPC_PEEK(this->contactMesh.vertices.size());

        for (UInt i = 0; i != patchCnt + 1; ++i)
        {
            BSIPC_WARN(
                "patch {}: vertex sep: {}, edge sep: {}, surface sep: {}",
                i,
                mesh.patchVertexSep[i],
                mesh.patchEdgeSep[i],
                mesh.patchSurfaceSep[i]
            );
        }

        this->contactMesh.friction = this->config.frictionMu;
        this->simlTime = 0;
    }

    void Solver::InitHessCache()
    {
        // Layout:
        //  Potential Energy Hessian | Inertia Hessian | Env Barrier Hessian | Seam Penalty Hessian

        UInt size = this->bsVertCntSum;

        UInt patchCnt = 0;
        for (const BSTargetInfo& info : this->targets)
            patchCnt += info.target->UDomainMax() * info.target->VDomainMax();

        // Hessian only consists of diagonal entries of the mass matrix. Mass matrix has size 3 * size == ControlPointCnt
        const UInt inertiaHessEntryCnt = 3 * size;

        // Flatten the local Hessians. Every hessian has 3-by-3 support, resulting in a 27-by-27 block
        // TODO check whether this size can be further reduced, by analyzing the structure of the Hessian
        const UInt potHessEntryCnt = 729 * patchCnt;
        //const UInt potHessEntryCnt = patchCnt * 15 * 2 * 9;

        this->curIndexInHessCache = 0;

        this->inertiaHessEntryCnt = inertiaHessEntryCnt;
        this->potHessEntryCnt = potHessEntryCnt;
        UInt membraneHessEntryCnt = potHessEntryCnt + inertiaHessEntryCnt;
        this->hessEntries.resize(membraneHessEntryCnt);
    }

    void Solver::PrecalculateBendingHess()
    {
        this->localBendingHess.resize(this->bsPatchCntSum);

        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            UInt bsPatchOffset = this->targets[idx].bsPatchIndexOffset.value();
            BSSurface* target = this->targets[idx].target;

            UInt patchCnt = target->GetPatchCnt();

            tbb::parallel_for(
                0, static_cast<Int>(patchCnt), 1, [&](Int idx)
                {
                    const QuadPoint& curBdQuadPoint = target->GetBdQuadPoints()[idx];
                    Float weight = curBdQuadPoint.Weight();
                    std::vector<UInt> localIndices = target->SupportedNodesAt(curBdQuadPoint);
                    Mat<27, 27> bdBlockHess = weight * Energy::BdHessBlockAt(curBdQuadPoint, localIndices, this->config, target);

                    this->localBendingHess[bsPatchOffset + idx] = bdBlockHess;
                }
            );
        }
    }

    void Solver::PrecomputeKKTComponents()
    {
        // System without animation sequence DoF
        UInt dofSize = 3 * this->bsVertCntSum;
        std::vector<std::unordered_map<UInt, Float>> seamSummedCoeffs;      // each seam point corresponds to an entry in the vector: vertex index -> coeff
        for (const SeamInfo& seamInfo : this->config.seamInfos)
        {
            BSSurface* surfI = this->targets[seamInfo.seamSurfIndices[0]].target;
            BSSurface* surfJ = this->targets[seamInfo.seamSurfIndices[1]].target;
            for (const SeamPoint& seamCoord : seamInfo.seamPoints)
            {
                Vec2 seamCoordI = seamCoord.first, seamCoordJ = seamCoord.second;

                std::vector<UInt> supportI = surfI->SupportedNodesAt(seamCoordI[0], seamCoordI[1]);
                std::vector<UInt> supportJ = surfJ->SupportedNodesAt(seamCoordJ[0], seamCoordJ[1]);

                UInt offsetI = this->targets[seamInfo.seamSurfIndices[0]].bsVertIndexOffset.value();
                UInt offsetJ = this->targets[seamInfo.seamSurfIndices[1]].bsVertIndexOffset.value();

                std::unordered_map<UInt, Float> curSeamPointCoeffs;

                // convention: coeff > 0 for surface I, coeff < 0 for surface J
                for (UInt i = 0; i != supportI.size(); ++i)
                {
                    UInt localIndex = supportI[i];
                    UInt globalIndex = offsetI + localIndex;
                    Float coeff = surfI->EvalCoeffAt(localIndex, seamCoordI[0], seamCoordI[1]);

                    if (std::abs(coeff) > 1e-4)
                    {
                        auto [it, inserted] = curSeamPointCoeffs.insert({ globalIndex, coeff });
                        if (!inserted)
                            it->second += coeff;
                    }
                }

                for (UInt j = 0; j != supportJ.size(); ++j)
                {
                    UInt localIndex = supportJ[j];
                    UInt globalIndex = offsetJ + localIndex;
                    Float coeff = -surfJ->EvalCoeffAt(localIndex, seamCoordJ[0], seamCoordJ[1]);

                    if (std::abs(coeff) > 1e-4)
                    {
                        auto [it, inserted] = curSeamPointCoeffs.insert({ globalIndex, coeff });
                        if (!inserted)
                            it->second += coeff;
                    }
                }

                seamSummedCoeffs.emplace_back(std::move(curSeamPointCoeffs));
            }
        }

        // Directly forming KKT system
        // Precomputed, using triplet for Hessian assmebly
        {
            SpMatData seamConstraintCoeffData;
            for (UInt i = 0; i != seamSummedCoeffs.size(); ++i)
            {
                const std::unordered_map<UInt, Float>& curSeamSummedCoeffs = seamSummedCoeffs[i];
                for (const auto& [idx, coeff] : curSeamSummedCoeffs)
                    for (UInt j = 0; j != 3; ++j)
                    {
                        seamConstraintCoeffData.emplace_back(3 * i + j, 3 * idx + j, coeff);
                    }
            }

            UInt hessSize = dofSize + 3 * seamSummedCoeffs.size();
            this->seamConstraintCnt = seamSummedCoeffs.size();

            this->seamConstraintCoeffs.resize(3 * seamSummedCoeffs.size(), dofSize);
            this->seamConstraintCoeffs.setFromTriplets(seamConstraintCoeffData.begin(), seamConstraintCoeffData.end());

            // See if need to precompute the extended matrix for animation sequence DoF
            if (this->animatedMeshCache.has_value())
            {
                UInt extendedSize = dofSize + 3 * this->animatedMeshCache.value().vertCnt;
                this->seamConstraintCoeffsExtended.resize(3 * seamSummedCoeffs.size(), extendedSize);
                this->seamConstraintCoeffsExtended.setFromTriplets(seamConstraintCoeffData.begin(), seamConstraintCoeffData.end());
            }
            else
                this->seamConstraintCoeffsExtended = this->seamConstraintCoeffs;
        }

        // Precompute the reordering and null-space projection matrix
        if (!seamSummedCoeffs.empty())
        {
            /// Precompute reordering matrix
            std::unordered_set<UInt> constrainedIndicesSet;
            for (const auto& seamCoeffs : seamSummedCoeffs)
                for (const auto& [idx, coeff] : seamCoeffs)
                    constrainedIndicesSet.insert(idx);

            UInt freeVertCnt = this->bsVertCntSum - static_cast<UInt>(constrainedIndicesSet.size());
            std::vector<UInt> indicesReorderMap;
            std::vector<UInt> constrainedIndicesVec(constrainedIndicesSet.begin(), constrainedIndicesSet.end());
            std::sort(constrainedIndicesVec.begin(), constrainedIndicesVec.end());
            indicesReorderMap.resize(this->bsVertCntSum);

            UInt reorderOffset = 0;
            for (UInt i = 0; i != this->bsVertCntSum; ++i)
            {
                if (constrainedIndicesSet.contains(i))
                    ++reorderOffset;
                else
                    indicesReorderMap[i - reorderOffset] = i;
            }
            for (UInt i = freeVertCnt; i != this->bsVertCntSum; ++i)
                indicesReorderMap[i] = constrainedIndicesVec[i - freeVertCnt];

            SpMatData rowReorderMatData;
            rowReorderMatData.reserve(3 * this->bsVertCntSum);

            for (UInt i = 0; i != this->bsVertCntSum; ++i)
                for (UInt j = 0; j != 3; ++j)
                    rowReorderMatData.emplace_back(3 * i + j, 3 * indicesReorderMap[i] + j, 1.);

            this->rowReorderingMat.resize(3 * this->bsVertCntSum, 3 * this->bsVertCntSum);
            this->rowReorderingMat.setFromTriplets(rowReorderMatData.begin(), rowReorderMatData.end());

            if (this->animatedMeshCache.has_value())
            {
                UInt animVertCnt = this->animatedMeshCache.value().vertCnt;
                UInt extendedDof = dofSize + 3 * animVertCnt;
                std::vector<UInt> indicesReorderMapExtended;
                indicesReorderMapExtended.resize(this->bsVertCntSum + animVertCnt);

                reorderOffset = 0;
                for (UInt i = 0; i != this->bsVertCntSum; ++i)
                {
                    if (constrainedIndicesSet.contains(i))
                        ++reorderOffset;
                    else
                        indicesReorderMapExtended[i - reorderOffset] = i;
                }
                for (UInt i = this->bsVertCntSum; i != this->bsVertCntSum + animVertCnt; ++i)
                    indicesReorderMapExtended[i - reorderOffset] = i;
                for (UInt i = freeVertCnt + animVertCnt; i != this->bsVertCntSum + animVertCnt; ++i)
                    indicesReorderMapExtended[i] = constrainedIndicesVec[i - freeVertCnt - animVertCnt];

                BSIPC_ASSERT(freeVertCnt == this->bsVertCntSum - reorderOffset, "[Reorder offset] and [free DoF Cnt] mismatch!");

                SpMatData rowReorderExtendedMatData;
                rowReorderExtendedMatData.reserve(3 * extendedDof);
                for (UInt i = 0; i != this->bsVertCntSum + animVertCnt; ++i)
                    for (UInt j = 0; j != 3; ++j)
                    {
                        rowReorderExtendedMatData.emplace_back(3 * i + j, 3 * indicesReorderMapExtended[i] + j, 1.);
                        if (3 * i + j >= 3 * extendedDof || 3 * indicesReorderMapExtended[i] + j >= 3 * extendedDof)
                            BSIPC_INFO("({}, {}) || {}", 3 * i + j, 3 * indicesReorderMapExtended[i] + j, 3 * extendedDof);
                    }

                this->rowReorderingMatExtended.resize(extendedDof, extendedDof);
                this->rowReorderingMatExtended.setFromTriplets(rowReorderExtendedMatData.begin(), rowReorderExtendedMatData.end());
            }
            else
                this->rowReorderingMatExtended = this->rowReorderingMat;

            /// Precompute nullspace projection matrix

            // Build C_c matrix on constrained DoF
            SpMatData constraintSelectioData;
            for (UInt i = 0; i != constrainedIndicesVec.size(); ++i)
            {
                UInt globalIdx = constrainedIndicesVec[i];
                for (UInt j = 0; j != 3; ++j)
                    constraintSelectioData.emplace_back(3 * globalIdx + j, 3 * i + j, 1.);
            }

            SpMat constraintSelectionMat;
            constraintSelectionMat.resize(dofSize, 3 * constrainedIndicesVec.size());
            constraintSelectionMat.setFromTriplets(constraintSelectioData.begin(), constraintSelectioData.end());

            SpMat condensedConstraintMat = this->seamConstraintCoeffs * constraintSelectionMat;                 // the Mat C_c

            //// Extract the projection matrix Z
            //Eigen::SparseQR<SpMat, Eigen::COLAMDOrdering<Int>> qrSolver;
            //qrSolver.compute(condensedConstraintMat.transpose());

            //if (qrSolver.info() != Eigen::Success)
            //    BSIPC_ERROR("QR decomposition on DoF reduction fails");

            //UInt qrRank = qrSolver.rank();
            //UInt constrainedDofSize = 3 * static_cast<UInt>(constrainedIndicesVec.size());
            //UInt nullspaceDim = constrainedDofSize - qrRank;
            //DMat qrQMat = qrSolver.matrixQ();
            //DMat Z = qrQMat.rightCols(nullspaceDim);

            //DMat prod = condensedConstraintMat * Z;
            //Float prodMax = prod.cwiseAbs().maxCoeff();
            //BSIPC_PEEK(prodMax);

            //BSIPC_INFO(
            //    "Partial DOF reduction: {} free + {} constrained -> {} free + {} reduced = {} total", 
            //    3 * freeVertCnt, constrainedDofSize, 
            //    3 * freeVertCnt, nullspaceDim, 
            //    3 * freeVertCnt + nullspaceDim
            //);

            //BSIPC_INFO("constrained Projection matrix [Z] dimension: {}x{}", Z.rows(), Z.cols());

            //// assemble the embedded projection matrix [I; Z]
            //SpMatData nullspaceProjMatData;
            //nullspaceProjMatData.reserve(3 * freeVertCnt + Z.rows() * Z.cols());

            //for (UInt i = 0; i != 3 * freeVertCnt; ++i)
            //    nullspaceProjMatData.emplace_back(i, i, 1.0);
            //for (UInt i = 0; i != Z.rows(); ++i)
            //    for (UInt j = 0; j != Z.cols(); ++j)
            //        if (std::abs(Z(i, j)) > 1e-8)
            //            nullspaceProjMatData.emplace_back(3 * freeVertCnt + i, 3 * freeVertCnt + j, Z(i, j));
            //BSIPC_WARN("3 * freeVertCnt: {}, Z: {}x{}, #nnz: {}", 3 * freeVertCnt, Z.rows(), Z.cols(), nullspaceProjMatData.size());

            //{
                // Extract the projection matrix Z using Gaussian elimination
                DMat C_dense = DMat(condensedConstraintMat);  // numConstraints x constrainedDofSize (no transpose)

                // Compute nullspace using FullPivLU
                Eigen::FullPivLU<DMat> lu(C_dense);
                UInt qrRank = lu.rank();
                UInt constrainedDofSize = 3 * static_cast<UInt>(constrainedIndicesVec.size());
                UInt nullspaceDim = constrainedDofSize - qrRank;

                if (nullspaceDim == 0)
                    BSIPC_ERROR("Nullspace dimension is zero - constraints are overconstrained");

                // Get kernel - this gives nullspace in the permuted coordinate system
                DMat Z_kernel = lu.kernel();  // size: (constrainedDofSize - rank) x nullspaceDim

                // Apply inverse column permutation to get back to original coordinates
                Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = lu.permutationQ();
                DMat Z_temp = DMat::Zero(constrainedDofSize, nullspaceDim);

                // The kernel corresponds to the free variables after permutation
                // Place kernel in the free variable positions (last nullspaceDim rows after permutation)
                Z_temp.bottomRows(nullspaceDim).setIdentity();
                Z_temp.topRows(qrRank) = -lu.matrixLU().topLeftCorner(qrRank, qrRank).triangularView<Eigen::Upper>().solve(
                    lu.matrixLU().topRightCorner(qrRank, nullspaceDim)
                );

                // Unpermute to get Z in original coordinates
                DMat Z = P * Z_temp;

                DMat prod = condensedConstraintMat * Z;
                Float prodMax = prod.cwiseAbs().maxCoeff();
                BSIPC_PEEK(prodMax);

                BSIPC_INFO(
                    "Partial DOF reduction: {} free + {} constrained -> {} free + {} reduced = {} total",
                    3 * freeVertCnt, constrainedDofSize,
                    3 * freeVertCnt, nullspaceDim,
                    3 * freeVertCnt + nullspaceDim
                );

                BSIPC_INFO("constrained Projection matrix [Z] dimension: {}x{}", Z.rows(), Z.cols());

                // assemble the embedded projection matrix [I; Z]
                SpMatData nullspaceProjMatData;
                nullspaceProjMatData.reserve(3 * freeVertCnt + Z.rows() * Z.cols());

                for (UInt i = 0; i != 3 * freeVertCnt; ++i)
                    nullspaceProjMatData.emplace_back(i, i, 1.0);
                for (UInt i = 0; i != Z.rows(); ++i)
                    for (UInt j = 0; j != Z.cols(); ++j)
                        if (std::abs(Z(i, j)) > 1e-8)
                            nullspaceProjMatData.emplace_back(3 * freeVertCnt + i, 3 * freeVertCnt + j, Z(i, j));
                BSIPC_WARN("3 * freeVertCnt: {}, Z: {}x{}, #nnz: {}", 3 * freeVertCnt, Z.rows(), Z.cols(), nullspaceProjMatData.size());
            //}

            this->constrainedDofProjMat.resize(3 * freeVertCnt + Z.rows(), 3 * freeVertCnt + Z.cols());
            this->constrainedDofProjMat.setFromTriplets(nullspaceProjMatData.begin(), nullspaceProjMatData.end());

            if (this->animatedMeshCache.has_value())
            {
                SpMatData nullspaceProjMatExtendedData;
                UInt animVertCnt = this->animatedMeshCache.value().vertCnt;
                UInt extendedSize = 3 * freeVertCnt + 3 * animVertCnt;
                nullspaceProjMatExtendedData.reserve(extendedSize + Z.rows() * Z.cols());

                for (UInt i = 0; i != extendedSize; ++i)
                    nullspaceProjMatExtendedData.emplace_back(i, i, 1.0);
                for (UInt i = 0; i != Z.rows(); ++i)
                    for (UInt j = 0; j != Z.cols(); ++j)
                        if (std::abs(Z(i, j)) > 1e-8)
                            nullspaceProjMatExtendedData.emplace_back(extendedSize + i, extendedSize + j, Z(i, j));

                this->constrainedDofProjMatExtended.resize(extendedSize + Z.rows(), extendedSize + Z.cols());
                this->constrainedDofProjMatExtended.setFromTriplets(
                    nullspaceProjMatExtendedData.begin(),
                    nullspaceProjMatExtendedData.end()
                );
            }
            else
                this->constrainedDofProjMatExtended = this->constrainedDofProjMat;

        }
    }

    void Solver::PreAssembleElasticityHessTemplate()
    {
        std::unordered_map<uint64_t, std::vector<ElasticityHessBlockSource>> hessBlockSourceMap;
        std::vector<std::vector<ElasticityGradBlockSource>> gradSourceMap(this->bsVertCntSum);

        auto CompressKey = IPC::compressPair;
        auto DecompressKey = IPC::decompressPair;

        /// Traverse all patches of all surfaces, record its link to entry indices in both gradient & Hessian
        for (const BSTargetInfo& info : this->targets)
        {
            UInt dofOffset = info.bsVertIndexOffset.value();
            UInt patchOffset = info.bsPatchIndexOffset.value();
            UInt uMax = info.target->UDomainMax();
            UInt vMax = info.target->VDomainMax();

            // TODO check whether inner iteration should be i or j
            for (UInt j = 0; j != vMax; ++j)
                for (UInt i = 0; i != uMax; ++i)
                {
                    UInt patchIdx = patchOffset + j * uMax + i;
                    Vec2 patchCenter = Vec2(i + 0.5, j + 0.5);
                    std::vector<UInt> localIndices = info.target->SupportedNodesAt(patchCenter[0], patchCenter[1]);
                    UInt localIndicesSize = localIndices.size();

                    std::vector<UInt> globalIndices(localIndicesSize);
                    for (UInt k = 0; k != localIndicesSize; ++k)
                        globalIndices[k] = dofOffset + localIndices[k];

                    for (UInt k = 0; k != globalIndices.size(); ++k)
                        for (UInt l = 0; l != globalIndices.size(); ++l)
                            hessBlockSourceMap[CompressKey(globalIndices[k], globalIndices[l])].emplace_back(patchIdx, k, l);

                    //for (UInt k = 0; k != globalIndices.size(); ++k)
                        //gradSourceMap[globalIndices[k]].emplace_back(patchIdx, k);
                }
        }

        /// Assemble the Hessian template.
        std::vector<std::vector<UInt>> blockRowColIndices;
        blockRowColIndices.resize(this->bsVertCntSum);

        // Pass 1: Categorize all the block indices over column indices
        for (const auto& block : hessBlockSourceMap)
        {
            std::array<UInt, 2> rowCol = DecompressKey(block.first);
            UInt row = rowCol[0], col = rowCol[1];
            blockRowColIndices[col].emplace_back(row);
        }

        // Pass 2: Sort within each column, and record the outer starts
        tbb::parallel_for(
            0, static_cast<Int>(blockRowColIndices.size()), 1, [&](Int col)
            {
                std::sort(blockRowColIndices[col].begin(), blockRowColIndices[col].end());
            }
        );
        UInt hessSize = 3 * this->bsVertCntSum;
        this->elasticityHessTemplate.resize(hessSize, hessSize);
        Int* hessOuterStarts = this->elasticityHessTemplate.outerIndexPtr();
        hessOuterStarts[0] = 0;
        for (UInt col = 0; col != this->bsVertCntSum; ++col)
        {
            UInt nnzInThisCol = 3 * static_cast<UInt>(blockRowColIndices[col].size());
            for (UInt i = 0; i != 3; ++i)
            {
                UInt index = 3 * col + i + 1;
                hessOuterStarts[index] = hessOuterStarts[index - 1] + nnzInThisCol;
            }
        }
        for (UInt col = 3 * this->bsVertCntSum; col != hessSize; ++col)
            hessOuterStarts[col] = hessOuterStarts[col - 1];

        // Pass 3: Compute all [innerIndices], and record the sources of all data entries
        UInt nnzCount = hessOuterStarts[hessSize - 1];
        this->elasticityHessTemplate.reserve(nnzCount);
        this->elasticityHessSources.resize(nnzCount);
        Int* hessInnerIndices = this->elasticityHessTemplate.innerIndexPtr();
        Float* hessValues = this->elasticityHessTemplate.valuePtr();
        std::memset(hessValues, 0, nnzCount * sizeof(Float));

        for (UInt blockColIdx = 0; blockColIdx != blockRowColIndices.size(); ++blockColIdx)
        {
            UInt blockCntInCurCol = static_cast<UInt>(blockRowColIndices[blockColIdx].size());
            for (UInt i = 0; i != blockCntInCurCol; ++i)
            {
                UInt blockRowIdx = blockRowColIndices[blockColIdx][i];
                for (UInt j = 0; j != 3; ++j)               // for each column of the block
                    for (UInt k = 0; k != 3; ++k)           // entry at (k, j)
                    {
                        UInt nnzIndex = hessOuterStarts[3 * blockColIdx + j] + 3 * i + k;
                        hessInnerIndices[nnzIndex] = 3 * blockRowIdx + k;

                        UInt compressedBlockKey = CompressKey(blockRowIdx, blockColIdx);
                        const std::vector<ElasticityHessBlockSource>& sources = hessBlockSourceMap[compressedBlockKey];
                        for (const auto& source : sources)
                        {
                            ElasticityHessEntrySource entrySource{
                                .globalPatchIdx = source.globalPatchIdx,
                                .rowIdx = 3 * source.blockRowIdx + k,
                                .colIdx = 3 * source.blockColIdx + j
                            };
                            this->elasticityHessSources[nnzIndex].emplace_back(entrySource);
                        }
                    }
            }
        }

        this->elasticityHessTemplate.resizeNonZeros(nnzCount);
        
        // Pass 4: Pre-fill the mass matrix, and [memset] all other entries to 0
        for (UInt i = 0; i != this->massMatDiagEntries.size(); ++i)
            for (UInt j = 0; j != 3; ++j)
            {
                UInt row = 3 * i + j;
                UInt col = 3 * i + j;
                Float value = this->massMatDiagEntries[i];
                this->elasticityHessTemplate.coeffRef(row, col) = value;
            }

        this->elasticityHessTemplate.makeCompressed();
        PrintSparseMatrixStructure(this->elasticityHessTemplate, "Elasticity Hessian Template");

    }

    void Solver::FilterSeamContacts()
    {
        BSIPC_WARN("#trig vertices: {}", this->contactMesh.vertices.size());
        BSIPC_WARN("#surfBackLinks entry: {}", this->contactMesh.surfBackLinks.size());

        for (const SeamInfo& seamInfo : this->config.seamInfos)
        {
            std::array<UInt, 2> surfIdx = seamInfo.seamSurfIndices;
            std::array<BSTargetInfo, 2> surfInfos = {
                this->targets[surfIdx[0]], 
                this->targets[surfIdx[1]]
            };
            std::array<UInt, 2> surfContactTrigStartIdx = {
                surfInfos[0].contactTrigVertIndexOffset.value(), surfInfos[1].contactTrigVertIndexOffset.value()
            };

            BSIPC_INFO("start: {} // {}", surfContactTrigStartIdx[0], surfContactTrigStartIdx[1]);

            std::array<UInt, 4> surfResolutions = {
                surfInfos[0].target->GetResolutionX(), surfInfos[0].target->GetResolutionY(),
                surfInfos[1].target->GetResolutionX(), surfInfos[1].target->GetResolutionY()
            };

            std::array<UInt, 4> surfTrigSizes = {
                surfInfos[0].target->GetLinearResolutionX(), surfInfos[0].target->GetLinearResolutionY(),
                surfInfos[1].target->GetLinearResolutionX(), surfInfos[1].target->GetLinearResolutionY()
            };

            BSIPC_WARN("trig sizes: {}x{} // {}x{}", surfTrigSizes[0], surfTrigSizes[1], surfTrigSizes[2], surfTrigSizes[3]);

            std::array<Float, 4> linearParamStep;
            for (UInt i = 0; i != 4; ++i)
                linearParamStep[i] = surfResolutions[i] / static_cast<Float>(surfTrigSizes[i] - 1);

            enum class BorderType
            {
                Lower, Upper, None
            };

            for (const SeamPoint& seamPoint : seamInfo.seamPoints)
            {
                // Identify bounds
                std::array<BorderType, 4> seamBorderTypes;
                std::array<Vec2, 2> seamCoords = { seamPoint.first, seamPoint.second };
                for (UInt i = 0; i != 2; ++i)
                {
                    for (UInt j = 0; j != 2; ++j)
                    {
                        UInt coordIdx = 2 * i + j;
                        if (std::abs(seamCoords[i][j]) < 1e-6)
                            seamBorderTypes[coordIdx] = BorderType::Lower;
                        else if (std::abs(seamCoords[i][j] - (surfResolutions[coordIdx] - 2)) < 1e-6)
                            seamBorderTypes[coordIdx] = BorderType::Upper;
                        else
                            seamBorderTypes[coordIdx] = BorderType::None;
                    }

                    BSIPC_ASSERT(
                        seamBorderTypes[2 * i] != BorderType::None || seamBorderTypes[2 * i + 1] != BorderType::None,
                        "At least one coordinate of a seam point on each surface should be on the border"
                    );
                }

                // Collect involved contact vertices
                std::array<std::vector<UInt>, 2> involvedContactVerts;

                for (UInt i = 0; i != 2; ++i)
                {
                    std::vector<UInt> idxX, idxY;
                    if (seamBorderTypes[0] == BorderType::Lower)
                        idxX = { 0 };
                    else if (seamBorderTypes[0] == BorderType::Upper)
                        idxX = { surfTrigSizes[2 * i] - 1 };
                    else
                    {
                        UInt lowerBound = static_cast<UInt>(std::floor(seamCoords[i][0] / linearParamStep[0]));
                        BSIPC_ASSERT(lowerBound < surfTrigSizes[2 * i] - 1, "Left bound index out of range");
                        idxX = { lowerBound, lowerBound + 1 };
                    }

                    if (seamBorderTypes[1] == BorderType::Lower)
                        idxY = { 0 };
                    else if (seamBorderTypes[1] == BorderType::Upper)
                        idxY = { surfTrigSizes[2 * i + 1] - 1 };
                    else
                    {
                        UInt lowerBound = static_cast<UInt>(std::floor(seamCoords[i][1] / linearParamStep[1]));
                        BSIPC_ASSERT(lowerBound < surfTrigSizes[2 * i + 1] - 1, "Lower bound index out of range");
                        idxY = { lowerBound, lowerBound + 1 };
                    }

                    for (UInt y : idxY)
                        for (UInt x : idxX)
                        {
                            involvedContactVerts[i].push_back(surfContactTrigStartIdx[i] + y * surfTrigSizes[i] + x);
                            BSIPC_INFO("{} // {} * {} // {}", surfContactTrigStartIdx[i], y, surfTrigSizes[i], x);
                        }
                }

                // Add primitives containing the involved contact vertices to the filter list
                std::array<std::unordered_set<UInt>, 2> involvedSurfVerts;
                std::array<std::unordered_set<UInt>, 2> involvedSurfEdges;
                std::array<std::unordered_set<UInt>, 2> involvedSurfTrigs;

                const auto& mesh = this->contactMesh;

                for (UInt i = 0; i != 2; ++i)
                    for (UInt vertIdx : involvedContactVerts[i])
                    {
                        BSIPC_PEEK(vertIdx);

                        const auto& surfVerts = mesh.surfBackLinks[vertIdx].surfVerts;
                        const auto& surfEdges = mesh.surfBackLinks[vertIdx].surfEdges;
                        const auto& surfTrigs = mesh.surfBackLinks[vertIdx].surfaces;

                        for (UInt vIdx : surfVerts)
                            involvedSurfVerts[i].insert(vIdx);

                        for (UInt eIdx : surfEdges)
                            involvedSurfEdges[i].insert(eIdx);

                        for (UInt tIdx : surfTrigs)
                            involvedSurfTrigs[i].insert(tIdx);
                    }

                // Write to the filter list in contact mesh

                // [point-trig]
                for (UInt vertIdx : involvedSurfVerts[0])
                    for (UInt trigIdx : involvedSurfTrigs[1])
                        this->contactMesh.pointTrigFilterSet.emplace(IPC::compressPair(vertIdx, trigIdx));

                for (UInt vertIdx : involvedSurfVerts[1])
                    for (UInt trigIdx : involvedSurfTrigs[0])
                        this->contactMesh.pointTrigFilterSet.emplace(IPC::compressPair(vertIdx, trigIdx));

                // [edge-edge]
                for (UInt edgeIdxI : involvedSurfEdges[0])
                    for (UInt edgeIdxJ : involvedSurfEdges[1])
                        this->contactMesh.edgeEdgeFilterSet.emplace(IPC::compressPair(edgeIdxI, edgeIdxJ));

                // [edge-trig]
                for (UInt edgeIdx : involvedSurfEdges[0])
                    for (UInt trigIdx : involvedSurfTrigs[1])
                        this->contactMesh.edgeTrigFilterSet.emplace(IPC::compressPair(edgeIdx, trigIdx));

                /// TEST PRINTS ///
                for (const uint64_t vertTrig : this->contactMesh.pointTrigFilterSet)
                {
                    auto [vIdx, tIdx] = IPC::decompressPair(vertTrig);
                    BSIPC_INFO("filter point-trig: {} // {}", vIdx, tIdx);
                }

                for (const uint64_t edgeEdge : this->contactMesh.edgeEdgeFilterSet)
                {
                    auto [eIdxI, eIdxJ] = IPC::decompressPair(edgeEdge);
                    BSIPC_INFO("filter edge-edge: {} // {}", eIdxI, eIdxJ);
                }

                for (const uint64_t edgeTrig : this->contactMesh.edgeTrigFilterSet)
                {
                    auto [eIdx, tIdx] = IPC::decompressPair(edgeTrig);
                    BSIPC_INFO("filter edge-trig: {} // {}", eIdx, tIdx);
                }
            }
        }
    }

    void Solver::NormalizeSeamCoordinates()
    {
        for (SeamInfo& info : this->config.seamInfos)
        {
            if (info.inCanonicalCoord)
            {
                BSIPC_INFO("info.seamSurfIndices: {}, {} // this->targets.size: {}", info.seamSurfIndices[0], info.seamSurfIndices[1], this->targets.size());
                BSIPC_ASSERT(info.seamSurfIndices[0] < this->targets.size() && info.seamSurfIndices[1] < this->targets.size(),
                    "seam surface indices should be within the number of recorded surfaces");

                const BSSurface* surfI = this->targets[info.seamSurfIndices[0]].target;
                const BSSurface* surfJ = this->targets[info.seamSurfIndices[1]].target;

                UInt surfIUMax = surfI->UDomainMax();
                UInt surfIVMax = surfI->VDomainMax();
                UInt surfJUMax = surfJ->UDomainMax();
                UInt surfJVMax = surfJ->VDomainMax();

                for (SeamPoint& point : info.seamPoints)
                {
                    point.first[0] *= surfIUMax;
                    point.first[1] *= surfIVMax;
                    point.second[0] *= surfJUMax;
                    point.second[1] *= surfJVMax;

                    BSIPC_INFO("seam: ({}, {}) -> ({}, {})",
                        point.first[0], point.first[1], point.second[0], point.second[1]
                    );
                }
            }
        }
    }

#if not defined BSIPC_DISABLE_RENDERER
    void Solver::InitRenderMesh(std::vector<std::pair<UInt, UInt>> bsTrigSizes)
    {
        this->renderMesh.Clear();

        std::vector<ATR::Vertex> vertices;

        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            BSTargetInfo& info = this->targets[idx];
            info.renderTrigVertIndexOffset = vertices.size();

            UInt uSize = bsTrigSizes[idx].first, vSize = bsTrigSizes[idx].second;

            // Culling mode is set to front = counter-clockwise
            // Iterate over the vertices
            for (UInt j = 0; j != vSize; ++j)
            {
                for (UInt i = 0; i != uSize; ++i)
                {
                    Float u = static_cast<Float>(i) / config.perPatchMesh;
                    Float v = static_cast<Float>(j) / config.perPatchMesh;

                    ATR::Vec3 pos = TransformVec3(info.target->EvalAt(u, v));
                    ATR::Vec3 normal = TransformVec3(info.target->NormalAt(u, v));

                    ATR::Vertex vertex(pos, normal);
                    vertices.push_back(vertex);
                }
            }
        }

        this->renderMesh.SetVertices(vertices);

        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            UInt uSize = bsTrigSizes[idx].first, vSize = bsTrigSizes[idx].second;
            UInt offset = this->targets[idx].renderTrigVertIndexOffset.value();

            // Iterate over the blocks, with (i, j) the lower-left corner vertex
            for (UInt j = 0; j != vSize - 1; ++j)
            {
                for (UInt i = 0; i != uSize - 1; ++i)
                {
                    UInt lowerLeft = offset + j * uSize + i;
                    UInt lowerRight = lowerLeft + 1;
                    UInt upperLeft = lowerLeft + uSize;
                    UInt upperRight = upperLeft + 1;

                    this->renderMesh.FormTriangle(lowerLeft, lowerRight, upperLeft);
                    this->renderMesh.FormTriangle(lowerRight, upperRight, upperLeft);
                }
            }
        }
    }
#endif

    void Solver::UpdateControlPoints(const std::vector<Vec3>& newPos)
    {
        // The order matters here!
        // Update underlying B-spline surface
        for (UInt i = 0; i != this->targets.size(); ++i)
        {
            BSTargetInfo& info = this->targets[i];
            UInt bsOffset = info.bsVertIndexOffset.value();
            std::vector<Vec3> curBSNewPos(newPos.begin() + bsOffset, newPos.begin() + bsOffset + info.target->GetControlPointCnt());

            info.target->UpdateControlPoints(curBSNewPos);
        }
        
        if (this->config.enableContact && !this->contactMesh.vertices.empty())
        {
            this->UpdateContactMesh();
        }

        if (!this->patchOrderMesh.empty())
        {
            this->UpdatePatchOrderMesh();
        }

        if (this->includeAnimatedMeshInSystem)
            this->UpdateAnimatedMesh(this->animatedMeshCache.value().curPos);
    }

    void Solver::UpdateAnimatedMesh(const std::vector<Vec3>& pos)
    {
        UInt vertCnt = this->animatedMeshCache.value().vertCnt;
        BSIPC_ASSERT(pos.size() == vertCnt, "Animated mesh position size mismatch");

        tbb::parallel_for(
            0, static_cast<Int>(this->animatedMeshCache.value().vertCnt), 1, [&](Int i)
            {
                this->contactMesh.vertices[this->animatedMeshCache.value().startIdx + i] = pos[i];
            }
        );
    }

    void Solver::UpdateMovingDBC(BSIPC_OUT std::vector<Vec3>& newPos, UInt curStep)
    {
        for (UInt i = 0; i != this->targets.size(); ++i)
        {
            for (MovingNode& node : this->config.surfacesInfo[i].movingNodes)
            {
                UInt index = node.index;
                if (index > newPos.size())
                {
                    BSIPC_WARN("Moving Node [{}] exceeding surface node count", index);
                    continue;
                }

                Vec3 vel = node.velocity;

                if (curStep < node.untilFrame)
                    for (UInt j = 0; j != 3; ++j)
                        newPos[index][j] += this->config.timestep * vel[j];
            }
        }
    }

    void Solver::UpdateTrigMesh(BSIPC_OUT IPC::mesh3D& mesh) const
    {
        Float perPatchMesh = config.perPatchMesh;

        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            const BSTargetInfo& info = this->targets[idx];

            //UInt uSize = static_cast<UInt>(info.target->UDomainMax() * perPatchMesh + 1);
            //UInt vSize = static_cast<UInt>(info.target->VDomainMax() * perPatchMesh + 1);

            UInt uSize = static_cast<UInt>(info.target->GetLinearResolutionX());
            UInt vSize = static_cast<UInt>(info.target->GetLinearResolutionY());
            UInt offset = info.contactTrigVertIndexOffset.value();

            // Update contact mesh
            for (UInt j = 0; j != vSize; ++j)
            {
                for (UInt i = 0; i != uSize; ++i)
                {
                    //Float u = static_cast<Float>(i) / perPatchMesh;
                    //Float v = static_cast<Float>(j) / perPatchMesh;

                    Float u = static_cast<Float>(i) / static_cast<Float>(uSize - 1) * info.target->UDomainMax();
                    Float v = static_cast<Float>(j) / static_cast<Float>(vSize - 1) * info.target->VDomainMax();

                    Vec3 pos = info.target->EvalAt(u, v);
                    UInt index = offset + j * uSize + i;
                    mesh.vertices[index] = pos;
                }
            }
        }
    }

    void Solver::UpdateContactMesh()
    {
        this->UpdateTrigMesh(this->contactMesh);

        IPC::mesh3D& mesh = this->contactMesh;

        double length = 0;
        for (const auto& edg : mesh.surfEdges)
        {
            length += (mesh.vertices[edg.first] -
                mesh.vertices[edg.second]).norm();
        }
        mesh.averageEdgeLength = length / (mesh.surfEdges.size());
    }

    void Solver::UpdatePatchOrderMesh()
    {
        Float perPatchMesh = config.perPatchMesh;
        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            const BSTargetInfo& info = this->targets[idx];
            TriangularMesh& mesh = this->patchOrderMesh[idx];

            UInt uSize = static_cast<UInt>(info.target->UDomainMax() * perPatchMesh + 1);
            UInt vSize = static_cast<UInt>(info.target->VDomainMax() * perPatchMesh + 1);

            // Update contact mesh
            for (UInt i = 0; i != uSize; ++i)
            {
                for (UInt j = 0; j != vSize; ++j)
                {
                    Float u = static_cast<Float>(i) / perPatchMesh;
                    Float v = static_cast<Float>(j) / perPatchMesh;

                    Vec3 pos = info.target->EvalAt(u, v);
                    Vec3 normal = info.target->NormalAt(u, v);
                    mesh.UpdateVertex(j * uSize + i, pos);
                    mesh.UpdateNormal(j * uSize + i, normal);
                }
            }
        }
    }

    void Solver::UpdateContactActiveSets()
    {
        this->spatialHash.build(this->contactMesh, this->contactMesh.averageEdgeLength);
        this->spatialHash.calculateActiveSet(this->contactMesh);
        this->ground.calculateActiveSet(this->contactMesh);
    }

    void Solver::InitPatchOrderMesh(const std::vector<std::pair<UInt, UInt>>& bsTrigSizes)
    {
        UInt surfCnt = this->targets.size();
        this->patchOrderMesh.resize(surfCnt);
        for (auto& mesh : this->patchOrderMesh)
            mesh.Clear();

        UInt trigStartOffset = 0;
        UInt bsStartIndex = 0;
        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            std::vector<Vertex> vertices;
            std::vector<Vec3> normals;
            BSTargetInfo& info = this->targets[idx];
            info.contactTrigVertIndexOffset = trigStartOffset;
            info.bsVertIndexOffset = bsStartIndex;

            UInt uSize = bsTrigSizes[idx].first, vSize = bsTrigSizes[idx].second;

            UInt curTargetTrigVertCnt = 0;
            for (UInt i = 0; i != uSize; ++i)
            {
                for (UInt j = 0; j != vSize; ++j)
                {
                    Float u = static_cast<Float>(i) / config.perPatchMesh;
                    Float v = static_cast<Float>(j) / config.perPatchMesh;
                    Vec3 pos = info.target->EvalAt(u, v);
                    Vec3 normal = info.target->NormalAt(u, v);
                    std::vector<SourceInfo> sources;

                    //BSIPC_INFO("Vertex: {}, {}, {}", pos[0], pos[1], pos[2]);
                    IVec2 supportStart = info.target->SupportStartOf(u, v);

                    Float coeff = 0;

                    for (UInt ctrlI = supportStart[0]; ctrlI < supportStart[0] + 3U && ctrlI < info.target->GetResolutionX(); ++ctrlI)
                    {
                        for (UInt ctrlJ = supportStart[1]; ctrlJ < supportStart[1] + 3U && ctrlJ < info.target->GetResolutionY(); ++ctrlJ)
                        {
                            Float coeffAtNode = info.target->EvalCoeffAt(ctrlI, ctrlJ, u, v);
                            if (coeffAtNode != 0)
                            {
                                UInt globalIndex = info.target->GetControlPointIndex(ctrlI, ctrlJ);
                                Float coeffAtIJ = info.target->EvalCoeffAt(ctrlI, ctrlJ, u, v);
                                coeff += coeffAtIJ;
                                sources.emplace_back(coeffAtIJ, globalIndex);
                                //BSIPC_INFO("- Source: {}, {}", coeffAtIJ, globalIndex);
                            }
                        }
                    }
                    ++curTargetTrigVertCnt;
                    ++trigStartOffset;

                    vertices.emplace_back(pos, sources);
                    normals.emplace_back(normal);
                }
            }
            info.trigVertCnt = curTargetTrigVertCnt;
            bsStartIndex += info.target->GetControlPointCnt();

            TriangularMesh& targetMesh = this->patchOrderMesh[idx];

            targetMesh.SetVertices(vertices);
            targetMesh.SetNormals(normals);

            for (UInt u = 0; u != info.target->UDomainMax(); ++u)
                for (UInt v = 0; v != info.target->VDomainMax(); ++v)
                    for (UInt i = 0; i != this->config.perPatchMesh; ++i)
                        for (UInt j = 0; j != this->config.perPatchMesh; ++j)
                        {
                            UInt lowerLeft = (v * this->config.perPatchMesh + j) * uSize + (u * this->config.perPatchMesh + i);
                            UInt lowerRight = lowerLeft + 1;
                            UInt upperLeft = lowerLeft + uSize;
                            UInt upperRight = upperLeft + 1;

                            UInt index = (u + v) * this->config.perPatchMesh + i + j;

                            if (index % 2)
                            {
                                targetMesh.FormTriangle(lowerLeft, lowerRight, upperLeft);
                                targetMesh.FormTriangle(lowerRight, upperRight, upperLeft);
                            }
                            else
                            {
                                targetMesh.FormTriangle(lowerLeft, lowerRight, upperRight);
                                targetMesh.FormTriangle(lowerLeft, upperRight, upperLeft);
                            }
                        }
        }
    }

#if not defined BSIPC_DISABLE_RENDERER
    void Solver::UpdateRenderMesh()
    {
        for (UInt idx = 0; idx != this->targets.size(); ++idx)
        {
            const BSTargetInfo& info = this->targets[idx];

            Float perPatchMesh = config.perPatchMesh;
            UInt uSize = static_cast<UInt>(info.target->UDomainMax() * perPatchMesh + 1);
            UInt vSize = static_cast<UInt>(info.target->VDomainMax() * perPatchMesh + 1);

            for (UInt j = 0; j != vSize; ++j)
            {
                for (UInt i = 0; i != uSize; ++i)
                {
                    Float u = static_cast<Float>(i) / perPatchMesh;
                    Float v = static_cast<Float>(j) / perPatchMesh;

                    ATR::Vec3 pos = TransformVec3(info.target->EvalAt(u, v));
                    ATR::Vec3 normal = TransformVec3(info.target->NormalAt(u, v));

                    UInt offset = info.renderTrigVertIndexOffset.value();

                    this->renderMesh.UpdateVertex(offset + j * uSize + i, pos, normal);
                }
            }
        }

        this->renderer->UpdateMesh(this->renderMesh);
    }
#endif
}
