#include "bspch.h"
#include "SeamPenalty.h"

namespace BSIPC
{
    namespace Energy
    {
        Float SeamValAt(const SeamPoint& seam, const SolverConfig& config, const BSSurface* surfI, const BSSurface* surfJ)
        {
            const Vec2& seamCoordI = seam.first;
            const Vec2& seamCoordJ = seam.second;

            Vec3 posI = surfI->EvalAt(seamCoordI[0], seamCoordI[1]);
            Vec3 posJ = surfJ->EvalAt(seamCoordJ[0], seamCoordJ[1]);

            Vec3 diff = posI - posJ;
            Float distance = diff.norm() * diff.norm();

            Float energy = config.seamStrength * distance;

            return energy;
        }

        Float SeamVal(const SolverConfig& config, const std::vector<BSTargetInfo>& infos)
        {
            const std::vector<SeamInfo> seams = config.seamInfos;

            std::vector<const BSSurface*> surfaces;
            for (const auto& target : infos)
                surfaces.push_back(target.target);

            Float seamEnergy = 0.0f;

            for (const auto& seam : seams)
            {
                BSIPC_ASSERT(seam.seamSurfIndices[0] < surfaces.size() && seam.seamSurfIndices[1] < surfaces.size(), 
                    "Invalid seam surface index");

                const BSSurface* surfI = surfaces[seam.seamSurfIndices[0]];
                const BSSurface* surfJ = surfaces[seam.seamSurfIndices[1]];

                for (const SeamPoint& seamCoord : seam.seamPoints)
                    seamEnergy += SeamValAt(seamCoord, config, surfI, surfJ);
            }

            return seamEnergy;
        }

        Mat<54, 1> SeamGradAt(const SeamPoint& seam, const SolverConfig& config,
            const BSSurface* surfI, const BSSurface* surfJ,
            std::vector<UInt> localIndicesI, std::vector<UInt> localIndicesJ)
        {
            Vec2 seamCoordI = seam.first;
            Vec2 seamCoordJ = seam.second;

            Vec3 posI = surfI->EvalAt(seamCoordI[0], seamCoordI[1]);
            Vec3 posJ = surfJ->EvalAt(seamCoordJ[0], seamCoordJ[1]);

            Vec3 diff = posI - posJ;

            Mat<54, 1> grad = Mat<54, 1>::Zero();

            for (UInt i = 0; i < localIndicesI.size(); ++i)
            {
                UInt localIndex = localIndicesI[i];
                Float localCoeff = surfI->EvalCoeffAt(localIndex, seamCoordI[0], seamCoordI[1]);
                Vec3 localGrad = 2 * config.seamStrength * diff * localCoeff;
                grad.block<3, 1>(i * 3, 0) += localGrad;
            }

            for (UInt j = 0; j < localIndicesJ.size(); ++j)
            {
                UInt localIndex = localIndicesJ[j];
                Float localCoeff = surfJ->EvalCoeffAt(localIndex, seamCoordJ[0], seamCoordJ[1]);
                Vec3 localGrad = -2 * config.seamStrength * diff * localCoeff;
                grad.block<3, 1>(j * 3 + 27, 0) += localGrad;
            }

            return grad;
        }

        void SeamGrad(std::vector<BSTargetInfo> targets, const SolverConfig& config, BSIPC_OUT DMat& seamGrad)
        {
            std::vector<BSTargetInfo> infos = targets;

            UInt size = 0;
            for (const auto& target : infos)
                size += target.target->GetControlPointCnt() * 3;

            seamGrad = DMat::Zero(size, 1);

            const std::vector<SeamInfo>& seamInfos = config.seamInfos;

            for (const SeamInfo& seamInfo : seamInfos)
            {
                BSIPC_ASSERT(seamInfo.seamSurfIndices[0] < infos.size() && seamInfo.seamSurfIndices[1] < infos.size(),
                    "Invalid seam surface index");

                BSSurface* surfI = infos[seamInfo.seamSurfIndices[0]].target;
                BSSurface* surfJ = infos[seamInfo.seamSurfIndices[1]].target;

                for (const SeamPoint& seamCoord : seamInfo.seamPoints)
                {
                    Vec2 seamCoordI = seamCoord.first;
                    Vec2 seamCoordJ = seamCoord.second;

                    std::vector<UInt> localIndicesI = surfI->SupportedNodesAt(seamCoordI[0], seamCoordI[1]);
                    std::vector<UInt> localIndicesJ = surfJ->SupportedNodesAt(seamCoordJ[0], seamCoordJ[1]);

                    UInt offsetI = infos[seamInfo.seamSurfIndices[0]].bsVertIndexOffset.value();
                    UInt offsetJ = infos[seamInfo.seamSurfIndices[1]].bsVertIndexOffset.value();

                    Mat<54, 1> localGrad = Energy::SeamGradAt(seamCoord, config, surfI, surfJ, localIndicesI, localIndicesJ);

#if defined BSIPC_NUMERIC_TEST
                    //DMat numericalSeamGrad = Energy::NumericalSeamGradAt(seamCoord, this->config, surfI, surfJ, localIndicesI, localIndicesJ);

                    //DMat seamGrad_T = localGrad.transpose();
                    //DMat numericalSeamGrad_T = numericalSeamGrad.transpose();
                    //DMat diff = seamGrad_T - numericalSeamGrad_T;

                    //this->AppendToLog("===== Seam Grad  Test =====\n\n");
                    //this->AppendToLog("Analytic Seam Grad: \n" + ToStr(seamGrad_T) + "\n");
                    //this->AppendToLog("Numerical Seam Grad: \n" + ToStr(numericalSeamGrad_T) + "\n");
                    //this->AppendToLog("Diff: \n" + ToStr(diff) + "\n");
                    //this->AppendToLog("\n<<<<< End Test >>>>>\n");
#endif
                    Energy::FillSeamGrad(localGrad, seamGrad, localIndicesI, localIndicesJ, offsetI, offsetJ);
                }
            }
        }

        void FillSeamGrad(const Mat<54, 1>& localGrad, BSIPC_OUT DMat& grad,
            std::vector<UInt> localIndicesI, std::vector<UInt> localIndicesJ, UInt offsetI, UInt offsetJ)
        {
            for (UInt i = 0; i != localIndicesI.size(); ++i)
            {
                UInt localIndex = localIndicesI[i];
                UInt globalIndex = localIndex + offsetI;
                grad.block<3, 1>(globalIndex * 3, 0) += localGrad.block<3, 1>(i * 3, 0);
            }

            for (UInt j = 0; j != localIndicesJ.size(); ++j)
            {
                UInt localIndex = localIndicesJ[j];
                UInt globalIndex = localIndex + offsetJ;
                grad.block<3, 1>(globalIndex * 3, 0) += localGrad.block<3, 1>(j * 3 + 27, 0);
            }
        }

        void SeamHessAt(const SeamPoint& seam, const SolverConfig& config,
            const BSSurface* surfI, const BSSurface* surfJ,
            std::vector<UInt> localIndicesI, std::vector<UInt> localIndicesJ,
            BSIPC_OUT DMat& localHess)
        {
            localHess = DMat::Zero(54, 54);

            Vec2 seamCoordI = seam.first, seamCoordJ = seam.second;

            for (UInt i = 0; i != localIndicesI.size(); ++i)
            {
                Float coeffI = surfI->EvalCoeffAt(localIndicesI[i], seamCoordI[0], seamCoordI[1]);

                for (UInt j = i; j != localIndicesI.size(); ++j)
                {
                    Float coeffJ = surfI->EvalCoeffAt(localIndicesI[j], seamCoordI[0], seamCoordI[1]);

                    for (UInt k = 0; k != 3; ++k)
                    {
                        localHess(i * 3 + k, j * 3 + k) = 2 * config.seamStrength * coeffI * coeffJ;
                        if (i != j)
                            localHess(j * 3 + k, i * 3 + k) = 2 * config.seamStrength * coeffI * coeffJ;
                    }
                }

                for (UInt j = 0; j != localIndicesJ.size(); ++j)
                {
                    Float coeffJ = surfJ->EvalCoeffAt(localIndicesJ[j], seamCoordJ[0], seamCoordJ[1]);

                    for (UInt k = 0; k != 3; ++k)
                    {
                        localHess(i * 3 + k, (9 + j) * 3 + k) = -2 * config.seamStrength * coeffI * coeffJ;
                        localHess((9 + j) * 3 + k, i * 3 + k) = -2 * config.seamStrength * coeffI * coeffJ;
                    }
                }
            }

            for (UInt i = 0; i != localIndicesJ.size(); ++i)
            {
                Float coeffI = surfJ->EvalCoeffAt(localIndicesJ[i], seamCoordJ[0], seamCoordJ[1]);

                for (UInt j = i; j != localIndicesJ.size(); ++j)
                {
                    Float coeffJ = surfJ->EvalCoeffAt(localIndicesJ[j], seamCoordJ[0], seamCoordJ[1]);

                    for (UInt k = 0; k != 3; ++k)
                    {
                        localHess((9 + i) * 3 + k, (9 + j) * 3 + k) = 2 * config.seamStrength * coeffI * coeffJ;
                        if (i != j)
                            localHess((9 + j) * 3 + k, (9 + i) * 3 + k) = 2 * config.seamStrength * coeffI * coeffJ;
                    }
                }
            }
            //localHess = SPDProjection(localHess);
        }

#if defined BSIPC_NUMERIC_TEST
        DMat NumericalSeamGradAt(const SeamPoint& seam, const SolverConfig& config, BSSurface* surfI, BSSurface* surfJ, 
            std::vector<UInt> localIndicesI, std::vector<UInt> localIndicesJ)
        {
            DMat result = DMat::Zero(54, 1);

            const Float h = 1e-6;

            for (UInt i = 0; i != localIndicesI.size(); ++i)
            {
                for (UInt coordI = 0; coordI != 3; ++coordI)
                {
                    UInt localIndex = localIndicesI[i];

                    Vec2 seamCoordI = seam.first;
                    Vec2 seamCoordJ = seam.second;

                    surfI->GetControlPoint(localIndex)[coordI] += h;
                    Float val1 = SeamValAt(seam, config, surfI, surfJ);

                    surfI->GetControlPoint(localIndex)[coordI] -= 2 * h;
                    Float val2 = SeamValAt(seam, config, surfI, surfJ);

                    surfI->GetControlPoint(localIndex)[coordI] += h;
                    result(3 * i + coordI, 0) = (val1 - val2) / (2 * h);
                }
            }

            for (UInt j = 0; j != localIndicesJ.size(); ++j)
            {
                for (UInt coordJ = 0; coordJ != 3; ++coordJ)
                {
                    UInt localIndex = localIndicesJ[j];

                    Vec2 seamCoordI = seam.first;
                    Vec2 seamCoordJ = seam.second;

                    surfJ->GetControlPoint(localIndex)[coordJ] += h;
                    Float val1 = SeamValAt(seam, config, surfI, surfJ);

                    surfJ->GetControlPoint(localIndex)[coordJ] -= 2 * h;
                    Float val2 = SeamValAt(seam, config, surfI, surfJ);

                    surfJ->GetControlPoint(localIndex)[coordJ] += h;
                    result(27 + 3 * j + coordJ, 0) = (val1 - val2) / (2 * h);
                }
            }

            return result;
        }

        DMat NumericalSeamHessAt(const SeamPoint& seam, const SolverConfig& config, BSSurface* surfI, BSSurface* surfJ,
            std::vector<UInt> localIndicesI, std::vector<UInt> localIndicesJ)
        {
            DMat hess = DMat::Zero(54, 54);

            const Float h = 1e-6;
            const Float hSq = h * h;

            std::vector<Vec3*> controlPtsHandle;
            for (UInt i = 0; i != localIndicesI.size(); ++i)
            {
                UInt localIndex = localIndicesI[i];
                controlPtsHandle.push_back(&surfI->GetControlPoint(localIndex));
            }
            for (UInt j = 0; j != localIndicesJ.size(); ++j)
            {
                UInt localIndex = localIndicesJ[j];
                controlPtsHandle.push_back(&surfJ->GetControlPoint(localIndex));
            }

            BSIPC_ASSERT(controlPtsHandle.size() == 18, "Wrong size in seam testing local Hessian");

            for (UInt i = 0; i != controlPtsHandle.size(); ++i)
                for (UInt j = 0; j != controlPtsHandle.size(); ++j)
                    for (UInt coordI = 0; coordI != 3; ++coordI)
                        for (UInt coordJ = 0; coordJ != 3; ++coordJ)
                        {
                            UInt indexI = 3 * i + coordI;
                            UInt indexJ = 3 * j + coordJ;

                            Vec3* controlPtI = controlPtsHandle[i];
                            Vec3* controlPtJ = controlPtsHandle[j];

                            (*controlPtI)[coordI] += h;
                            (*controlPtJ)[coordJ] += h;
                            Float val1 = SeamValAt(seam, config, surfI, surfJ);

                            (*controlPtI)[coordI] -= 2 * h;
                            Float val2 = SeamValAt(seam, config, surfI, surfJ);

                            (*controlPtI)[coordI] += 2 * h;
                            (*controlPtJ)[coordJ] -= 2 * h;
                            Float val3 = SeamValAt(seam, config, surfI, surfJ);

                            (*controlPtI)[coordI] -= 2 * h;
                            Float val4 = SeamValAt(seam, config, surfI, surfJ);

                            (*controlPtI)[coordI] += h;
                            (*controlPtJ)[coordJ] += h;

                            hess(indexI, indexJ) = (val1 - val2 - val3 + val4) / (4 * hSq);
                        }

            return hess;
        }

        DMat NumericalSeamGrad(const SolverConfig& config, const std::vector<BSTargetInfo>& infos)
        {
            const Float h = 1e-6;

            UInt bsNodeCnt = 0;
            for (const BSTargetInfo& target : infos)
                bsNodeCnt += target.target->GetControlPointCnt();

            DMat result = DMat::Zero(bsNodeCnt * 3, 1);

            for (UInt idx = 0; idx != infos.size(); ++idx)
            {
                UInt bsOffset = infos[idx].bsVertIndexOffset.value();
                BSSurface* target = infos[idx].target;

                for (UInt i = 0; i != target->GetControlPointCnt(); ++i)
                {
                    Vec3& controlPt = target->GetControlPoint(i);

                    for (UInt coord = 0; coord != 3; ++coord)
                    {
                        controlPt[coord] += h;
                        Float val1 = SeamVal(config, infos);

                        controlPt[coord] -= 2 * h;
                        Float val2 = SeamVal(config, infos);

                        controlPt[coord] += h;
                        result((bsOffset + i) * 3 + coord, 0) = (val1 - val2) / (2 * h);
                    }
                }
            }

            return result;
        }
        

        DMat NumericalSeamHess(const SolverConfig& config, const std::vector<BSTargetInfo>& infos)
        {
            const Float h = 1e-6;
            const Float hSq = h * h;

            UInt bsNodeCnt = 0;
            for (const BSTargetInfo& target : infos)
                bsNodeCnt += target.target->GetControlPointCnt();

            DMat hess = DMat::Zero(bsNodeCnt * 3, bsNodeCnt * 3);

            for (UInt idxI = 0; idxI != infos.size(); ++idxI)
            {
                UInt bsOffsetI = infos[idxI].bsVertIndexOffset.value();
                BSSurface* targetI = infos[idxI].target;
                for (UInt idxJ = 0; idxJ != infos.size(); ++idxJ)
                {
                    UInt bsOffsetJ = infos[idxJ].bsVertIndexOffset.value();
                    BSSurface* targetJ = infos[idxJ].target;
                    for (UInt i = 0; i != targetI->GetControlPointCnt(); ++i)
                    {
                        Vec3& controlPtI = targetI->GetControlPoint(i);
                        for (UInt j = 0; j != targetJ->GetControlPointCnt(); ++j)
                        {
                            Vec3& controlPtJ = targetJ->GetControlPoint(j);
                            for (UInt coordI = 0; coordI != 3; ++coordI)
                                for (UInt coordJ = 0; coordJ != 3; ++coordJ)
                                {
                                    UInt indexI = (bsOffsetI + i) * 3 + coordI;
                                    UInt indexJ = (bsOffsetJ + j) * 3 + coordJ;

                                    controlPtI[coordI] += h;
                                    controlPtJ[coordJ] += h;
                                    Float val1 = SeamVal(config, infos);

                                    controlPtI[coordI] -= 2 * h;
                                    Float val2 = SeamVal(config, infos);

                                    controlPtI[coordI] += 2 * h;
                                    controlPtJ[coordJ] -= 2 * h;
                                    Float val3 = SeamVal(config, infos);

                                    controlPtI[coordI] -= 2 * h;
                                    Float val4 = SeamVal(config, infos);

                                    controlPtI[coordI] += h;
                                    controlPtJ[coordJ] += h;

                                    hess(indexI, indexJ) = (val1 - val2 - val3 + val4) / (4 * hSq);
                                }
                        }
                    }
                }
            }

            return hess;
        }
#endif

    }
}
