#include "bspch.h"

#include "Bending.h"

#include "Utility/EigenHelper.h"

namespace BSIPC
{
    namespace Energy
    {
        Float BdValAt(const QuadPoint& quadPoint, const SolverConfig& config, const BSSurface* const surface)
        {
            Float u = quadPoint.U(), v = quadPoint.V();

            BSIPC_ASSERT(u >= 0 && u <= surface->UDomainMax(), "u out of bounds");
            BSIPC_ASSERT(v >= 0 && v <= surface->VDomainMax(), "v out of bounds");

            Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(surface->GetResolutionX() - 3));
            Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(surface->GetResolutionY() - 3));
            std::array<UInt, 3> uOffset = surface->SelectSupportedOffset(u);
            std::array<UInt, 3> vOffset = surface->SelectSupportedOffset(v);

            std::array<BSBasisFn, 3> horzValFn = surface->SelectSupportedBasisFn(u, surface->GetResolutionX(), BSBasis::BasisFns);
            std::array<BSBasisFn, 3> vertValFn = surface->SelectSupportedBasisFn(v, surface->GetResolutionY(), BSBasis::BasisFns);

            std::array<BSBasisFn, 3> horzDerivFn = surface->SelectSupportedBasisFn(u, surface->GetResolutionX(), BSBasisDerivative::BasisFns);
            std::array<BSBasisFn, 3> vertDerivFn = surface->SelectSupportedBasisFn(v, surface->GetResolutionY(), BSBasisDerivative::BasisFns);

            std::array<BSBasisFn, 3> horzSecondDerivFn = surface->SelectSupportedBasisFn(u, surface->GetResolutionX(), BSBasisSecondDerivative::BasisFns);
            std::array<BSBasisFn, 3> vertSecondDerivFn = surface->SelectSupportedBasisFn(v, surface->GetResolutionY(), BSBasisSecondDerivative::BasisFns);

            // Laplacian at quadPoint (\Delta \phi(u, v))
            Vec3 LP_Pos = Vec3(0, 0, 0);

            // Cached PDs
            const Mat<2, 2>& pd_P_M = quadPoint.PD_Param_Material();
            const Mat<2, 3>& pd2_P_M = quadPoint.PD2_Param_Material();

            Float ux = pd_P_M(0, 0), uy = pd_P_M(0, 1), vx = pd_P_M(1, 0), vy = pd_P_M(1, 1);
            Float uxx = pd2_P_M(0, 0), uyy = pd2_P_M(0, 2), vxx = pd2_P_M(0, 2), vyy = pd2_P_M(1, 2);

            for (UInt i = 0; i != 3; ++i)       // Support is 3-by-3
            {
                for (UInt j = 0; j != 3; ++j)
                {
                    Float uEval = (u - uFloor) + uOffset[i];
                    Float vEval = (v - vFloor) + vOffset[j];

                    Float horzVal = horzValFn[i](uEval);
                    Float vertVal = vertValFn[j](vEval);
                    Float horzDerivVal = horzDerivFn[i](uEval);
                    Float vertDerivVal = vertDerivFn[j](vEval);
                    Float horzSecondDerivVal = horzSecondDerivFn[i](uEval);
                    Float vertSecondDerivVal = vertSecondDerivFn[j](vEval);

                    Float coeff =
                        (ux * ux + uy * uy) * horzSecondDerivVal * vertVal +
                        (vx * vx + vy * vy) * horzVal * vertSecondDerivVal +
                        2 * (ux * vx + uy * vy) * horzDerivVal * vertDerivVal +
                        (uxx + uyy) * horzDerivVal * vertVal +
                        (vxx + vyy) * horzVal * vertDerivVal;

                    Vec3 curSumSegment = coeff * surface->GetControlPoint(uFloor + i, vFloor + j);
                    LP_Pos += curSumSegment;
                }
            }

            Float meanCurvatureSq = LP_Pos.dot(LP_Pos);
            Float result = 0.5 * config.bendingStiffness * meanCurvatureSq;

            //BSIPC_INFO("Bending energy at ({}, {}): {}", u, v, result);

            return result;
        }

        DMat BdGradBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface)
        {
            Float u = quadPoint.U(), v = quadPoint.V();

            // Cached PDs
            const Mat<2, 2>& pd_P_M = quadPoint.PD_Param_Material();
            const Mat<2, 3>& pd2_P_M = quadPoint.PD2_Param_Material();

            Float ux = pd_P_M(0, 0), uy = pd_P_M(0, 1), vx = pd_P_M(1, 0), vy = pd_P_M(1, 1);
            Float uxx = pd2_P_M(0, 0), uyy = pd2_P_M(0, 2), vxx = pd2_P_M(0, 2), vyy = pd2_P_M(1, 2);

            DMat bdGrad = DMat::Zero(27, 1);

            BSIPC_ASSERT(u >= 0 && u <= surface->UDomainMax(), "u out of bounds");
            BSIPC_ASSERT(v >= 0 && v <= surface->VDomainMax(), "v out of bounds");
            BSIPC_ASSERT(indices.size() == 9, "support must has 9 nodes");

            for (UInt baseIndex = 0; baseIndex != indices.size(); ++baseIndex)
            {
                IVec2 basePosIndex = surface->GetControlPointPosIndex(indices[baseIndex]);

                Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(surface->GetResolutionX() - 3));
                Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(surface->GetResolutionY() - 3));

                // Compute constant term c_ij
                Vec2 uSupportAtBase = surface->SingleDimSupport(basePosIndex[0], surface->GetResolutionX());
                Vec2 vSupportAtBase = surface->SingleDimSupport(basePosIndex[1], surface->GetResolutionY());

                BSBasisFn baseHorzValFn = surface->SelectBasisFnByIndex(basePosIndex[0], BSBasis::BasisFns, surface->GetResolutionX());
                BSBasisFn baseVertValFn = surface->SelectBasisFnByIndex(basePosIndex[1], BSBasis::BasisFns, surface->GetResolutionY());
                BSBasisFn baseHorzDerivFn = surface->SelectBasisFnByIndex(basePosIndex[0], BSBasisDerivative::BasisFns, surface->GetResolutionX());
                BSBasisFn baseVertDerivFn = surface->SelectBasisFnByIndex(basePosIndex[1], BSBasisDerivative::BasisFns, surface->GetResolutionY());
                BSBasisFn baseHorzDeriv2Fn = surface->SelectBasisFnByIndex(basePosIndex[0], BSBasisSecondDerivative::BasisFns, surface->GetResolutionX());
                BSBasisFn baseVertDeriv2Fn = surface->SelectBasisFnByIndex(basePosIndex[1], BSBasisSecondDerivative::BasisFns, surface->GetResolutionY());

                Float baseUEvalOffset = u - uSupportAtBase[0];
                Float baseVEvalOffset = v - vSupportAtBase[0];

                Float baseHorzVal = baseHorzValFn(baseUEvalOffset);
                Float baseVertVal = baseVertValFn(baseVEvalOffset);
                Float baseHorzDerivVal = baseHorzDerivFn(baseUEvalOffset);
                Float baseVertDerivVal = baseVertDerivFn(baseVEvalOffset);
                Float baseHorzDeriv2Val = baseHorzDeriv2Fn(baseUEvalOffset);
                Float baseVertDeriv2Val = baseVertDeriv2Fn(baseVEvalOffset);

                Float constCoeff =
                    (ux * ux + uy * uy) * baseHorzDeriv2Val * baseVertVal +
                    (vx * vx + vy * vy) * baseHorzVal * baseVertDeriv2Val +
                    2 * (ux * vx + uy * vy) * baseHorzDerivVal * baseVertDerivVal +
                    (uxx + uyy) * baseHorzDerivVal * baseVertVal +
                    (vxx + vyy) * baseHorzVal * baseVertDerivVal;

                // Compute summation of other terms \sum_{k, l} c_kl * P_kl
                std::array<UInt, 3> uOffset = surface->SelectSupportedOffset(u);
                std::array<UInt, 3> vOffset = surface->SelectSupportedOffset(v);

                std::array<BSBasisFn, 3> horzEval = surface->SelectSupportedBasisFn(u, surface->GetResolutionX(), BSBasis::BasisFns);
                std::array<BSBasisFn, 3> vertEval = surface->SelectSupportedBasisFn(v, surface->GetResolutionY(), BSBasis::BasisFns);
                std::array<BSBasisFn, 3> horzDeriv = surface->SelectSupportedBasisFn(u, surface->GetResolutionX(), BSBasisDerivative::BasisFns);
                std::array<BSBasisFn, 3> vertDeriv = surface->SelectSupportedBasisFn(v, surface->GetResolutionY(), BSBasisDerivative::BasisFns);
                std::array<BSBasisFn, 3> horzDeriv2 = surface->SelectSupportedBasisFn(u, surface->GetResolutionX(), BSBasisSecondDerivative::BasisFns);
                std::array<BSBasisFn, 3> vertDeriv2 = surface->SelectSupportedBasisFn(v, surface->GetResolutionY(), BSBasisSecondDerivative::BasisFns);

                // [ dH^2/dP_ijx dH^2/dP_ijy dH^2/dP_ijz ]
                // where (i, j) = basePosIndex
                Vec3 result = Vec3(0, 0, 0);
                for (UInt i = 0; i != 3; ++i)
                {
                    for (UInt j = 0; j != 3; ++j)
                    {
                        Float uEval = (u - uFloor) + uOffset[i];
                        Float vEval = (v - vFloor) + vOffset[j];

                        Float curSegmentCoeff =
                            (ux * ux + uy * uy) * horzDeriv2[i](uEval) * vertEval[j](vEval) +
                            (vx * vx + vy * vy) * horzEval[i](uEval) * vertDeriv2[j](vEval) +
                            2 * (ux * vx + uy * vy) * horzDeriv[i](uEval) * vertDeriv[j](vEval) +
                            (uxx + uyy) * horzDeriv[i](uEval) * vertEval[j](vEval) +
                            (vxx + vyy) * horzEval[i](uEval) * vertDeriv[j](vEval);

                        Vec3 step = curSegmentCoeff * surface->GetControlPoint(uFloor + i, vFloor + j);

                        IVec2 supportStart = surface->SupportStartOf(u, v);

                        result += 2. * constCoeff * step;
                    }
                }

                bdGrad.block<3, 1>(3 * baseIndex, 0) = result;
            }

            DMat scaledBdGrad = 0.5 * config.bendingStiffness * bdGrad;

            return scaledBdGrad;
        }

        Mat<27, 27> BdHessBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface)
        {
            Float u = quadPoint.U(), v = quadPoint.V();

            Mat<27, 27> bdHess = Mat<27, 27>::Zero();

            // Cached PDs
            const Mat<2, 2>& pd_P_M = quadPoint.PD_Param_Material();
            const Mat<2, 3>& pd2_P_M = quadPoint.PD2_Param_Material();

            Float ux = pd_P_M(0, 0), uy = pd_P_M(0, 1), vx = pd_P_M(1, 0), vy = pd_P_M(1, 1);
            Float uxx = pd2_P_M(0, 0), uyy = pd2_P_M(0, 2), vxx = pd2_P_M(0, 2), vyy = pd2_P_M(1, 2);

            BSIPC_ASSERT(u >= 0 && u <= surface->UDomainMax(), "u out of bounds");
            BSIPC_ASSERT(v >= 0 && v <= surface->VDomainMax(), "v out of bounds");
            BSIPC_ASSERT(indices.size() == 9, "support must has 9 nodes");

            for (UInt i = 0; i != 9; ++i)
            {
                for (UInt j = 0; j != 9; ++j)
                {
                    IVec2 node1 = surface->GetControlPointPosIndex(indices[i]);
                    IVec2 node2 = surface->GetControlPointPosIndex(indices[j]);

                    BSBasisFn horzEval1 = surface->SelectBasisFnByIndex(node1[0], BSBasis::BasisFns, surface->GetResolutionX());
                    BSBasisFn vertEval1 = surface->SelectBasisFnByIndex(node1[1], BSBasis::BasisFns, surface->GetResolutionY());
                    BSBasisFn horzEval2 = surface->SelectBasisFnByIndex(node2[0], BSBasis::BasisFns, surface->GetResolutionX());
                    BSBasisFn vertEval2 = surface->SelectBasisFnByIndex(node2[1], BSBasis::BasisFns, surface->GetResolutionY());

                    BSBasisFn horzDerivFn1 = surface->SelectBasisFnByIndex(node1[0], BSBasisDerivative::BasisFns, surface->GetResolutionX());
                    BSBasisFn vertDerivFn1 = surface->SelectBasisFnByIndex(node1[1], BSBasisDerivative::BasisFns, surface->GetResolutionY());
                    BSBasisFn horzDerivFn2 = surface->SelectBasisFnByIndex(node2[0], BSBasisDerivative::BasisFns, surface->GetResolutionX());
                    BSBasisFn vertDerivFn2 = surface->SelectBasisFnByIndex(node2[1], BSBasisDerivative::BasisFns, surface->GetResolutionY());

                    BSBasisFn horzDeriv2Fn1 = surface->SelectBasisFnByIndex(node1[0], BSBasisSecondDerivative::BasisFns, surface->GetResolutionX());
                    BSBasisFn vertDeriv2Fn1 = surface->SelectBasisFnByIndex(node1[1], BSBasisSecondDerivative::BasisFns, surface->GetResolutionY());
                    BSBasisFn horzDeriv2Fn2 = surface->SelectBasisFnByIndex(node2[0], BSBasisSecondDerivative::BasisFns, surface->GetResolutionX());
                    BSBasisFn vertDeriv2Fn2 = surface->SelectBasisFnByIndex(node2[1], BSBasisSecondDerivative::BasisFns, surface->GetResolutionY());

                    Vec2 uSupport1 = surface->SingleDimSupport(node1[0], surface->GetResolutionX());
                    Vec2 vSupport1 = surface->SingleDimSupport(node1[1], surface->GetResolutionY());
                    Vec2 uSupport2 = surface->SingleDimSupport(node2[0], surface->GetResolutionX());
                    Vec2 vSupport2 = surface->SingleDimSupport(node2[1], surface->GetResolutionY());

                    Float uEvalOffset1 = u - uSupport1[0];
                    Float vEvalOffset1 = v - vSupport1[0];
                    Float uEvalOffset2 = u - uSupport2[0];
                    Float vEvalOffset2 = v - vSupport2[0];

                    Float coeff1 =
                        (ux * ux + uy * uy) * horzDeriv2Fn1(uEvalOffset1) * vertEval1(vEvalOffset1) +
                        (vx * vx + vy * vy) * horzEval1(uEvalOffset1) * vertDeriv2Fn1(vEvalOffset1) +
                        2 * (ux * vx + uy * vy) * horzDerivFn1(uEvalOffset1) * vertDerivFn1(vEvalOffset1) +
                        (uxx + uyy) * horzDerivFn1(uEvalOffset1) * vertEval1(vEvalOffset1) +
                        (vxx + vyy) * horzEval1(uEvalOffset1) * vertDerivFn1(vEvalOffset1);

                    Float coeff2 =
                        (ux * ux + uy * uy) * horzDeriv2Fn2(uEvalOffset2) * vertEval2(vEvalOffset2) +
                        (vx * vx + vy * vy) * horzEval2(uEvalOffset2) * vertDeriv2Fn2(vEvalOffset2) +
                        2 * (ux * vx + uy * vy) * horzDerivFn2(uEvalOffset2) * vertDerivFn2(vEvalOffset2) +
                        (uxx + uyy) * horzDerivFn2(uEvalOffset2) * vertEval2(vEvalOffset2) +
                        (vxx + vyy) * horzEval2(uEvalOffset2) * vertDerivFn2(vEvalOffset2);

                    Float result = coeff1 * coeff2;

                    for (UInt ind = 0; ind != 3; ++ind)
                    {
                        UInt baseIndex1 = 3 * i + ind;
                        UInt baseIndex2 = 3 * j + ind;

                        bdHess(baseIndex1, baseIndex2) = 2. * result;
                    }
                }
            }

            return 0.5 * config.bendingStiffness * bdHess;
        }

#if defined BSIPC_NUMERIC_TEST
        DMat NumericalBdGradBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, BSIPC_OUT BSSurface* const target)
        {
            DMat grad = DMat::Zero(3 * indices.size(), 1);
            const Float h = 1e-8;

            for (UInt i = 0; i != indices.size(); ++i)
            {
                for (UInt index = 0; index != 3; ++index)
                {
                    UInt nodeIndex = indices[i];
                    Vec3& perturbVertex = target->GetControlPoint(nodeIndex);

                    perturbVertex[index] += h;
                    target->UpdateQuadPointCache();
                    Float upVal = BdValAt(quadPoint, config, target);

                    perturbVertex[index] -= 2 * h;
                    target->UpdateQuadPointCache();
                    Float downVal = BdValAt(quadPoint, config, target);

                    perturbVertex[index] += h;

                    grad(3 * i + index, 0) = (upVal - downVal) / (2 * h);
                }
            }

            return grad;
        }
#endif

#if defined BSIPC_NUMERIC_TEST 
        DMat NumericalBdHessBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, BSIPC_OUT BSSurface* const target)
        {
            const Float h = 1e-8;
            const Float hSq = h * h;
            DMat hess = DMat::Zero(27, 27);

            for (UInt i = 0; i != indices.size(); ++i)
                for (UInt j = 0; j != indices.size(); ++j)
                    for (UInt indexI = 0; indexI != 3; ++indexI)
                        for (UInt indexJ = 0; indexJ != 3; ++indexJ)
                        {
                            UInt matIndexI = 3 * i + indexI;
                            UInt matIndexJ = 3 * j + indexJ;

                            Vec3& perturbVertexI = target->GetControlPoint(indices[i]);
                            Vec3& perturbVertexJ = target->GetControlPoint(indices[j]);

                            // Up, f(... x_i + h ... x_j + h ...)
                            perturbVertexI[indexI] += h;
                            perturbVertexJ[indexJ] += h;
                            target->UpdateQuadPointCache();
                            Float upSsVal = BdValAt(quadPoint, config, target);

                            // Left, f(... x_i - h ... x_j + h ...)
                            perturbVertexI[indexI] -= 2 * h;
                            target->UpdateQuadPointCache();
                            Float leftSsVal = BdValAt(quadPoint, config, target);

                            // Right, f(... x_i + h ... x_j - h ...)
                            perturbVertexI[indexI] += 2 * h;
                            perturbVertexJ[indexJ] -= 2 * h;
                            target->UpdateQuadPointCache();
                            Float rightSsVal = BdValAt(quadPoint, config, target);

                            // Down, f(... x_i - h ... x_j - h ...)
                            perturbVertexI[indexI] -= 2 * h;
                            target->UpdateQuadPointCache();
                            Float downSsVal = BdValAt(quadPoint, config, target);

                            Float hessEntry = (upSsVal - rightSsVal - leftSsVal + downSsVal) / (4 * hSq);
                            hess(matIndexI, matIndexJ) = hessEntry;
                        }

            return hess;
        }

        DMat NumericalBdHess(const SolverConfig& config, const std::vector<BSTargetInfo>& infos)
        {
            UInt size = 0;

            for (const BSTargetInfo& info : infos)
                size += info.target->GetControlPointCnt();

            auto BdVal = [&]() -> Float
                {
                    Float bdEnergyVal = 0;

                    for (const BSTargetInfo& info : infos)
                    {
                        const std::vector<QuadPoint> quadPoints = info.target->GetBdQuadPoints();
                        const UInt iterSize = quadPoints.size();

                        Float curBdEnergyVal =
                            oneapi::tbb::parallel_reduce(
                                oneapi::tbb::blocked_range<UInt>(0, iterSize),
                                0.,
                                [&](const oneapi::tbb::blocked_range<UInt>& range, Float localSum) -> Float
                                {
                                    for (UInt i = range.begin(); i != range.end(); ++i)
                                    {
                                        const QuadPoint& curQuadPoint = quadPoints[i];

                                        Float weight = curQuadPoint.Weight();
                                        localSum += weight * BdValAt(curQuadPoint, config, info.target);
                                    }

                                    return localSum;
                                },
                                [&](Float a, Float b) -> Float
                                {
                                    return a + b;
                                }
                            );

                        bdEnergyVal += curBdEnergyVal;
                    }

                    return bdEnergyVal;
                };

            DMat hess = DMat::Zero(3 * size, 3 * size);

            const Float h = 1e-8;
            const Float hSq = h * h;

            for (const BSTargetInfo& info : infos)
                for (UInt i = 0; i != size; ++i)
                    for (UInt j = 0; j != size; ++j)
                        for (UInt indexI = 0; indexI != 3; ++indexI)
                            for (UInt indexJ = 0; indexJ != 3; ++indexJ)
                            {
                                UInt matIndexI = 3 * i + indexI;
                                UInt matIndexJ = 3 * j + indexJ;

                                Vec3& perturbVertexI = info.target->GetControlPoint(i);
                                Vec3& perturbVertexJ = info.target->GetControlPoint(j);

                                // Up, f(... x_i + h ... x_j + h ...)
                                perturbVertexI[indexI] += h;
                                perturbVertexJ[indexJ] += h;
                                info.target->UpdateQuadPointCache();
                                Float upSsVal = BdVal();

                                // Left, f(... x_i - h ... x_j + h ...)
                                perturbVertexI[indexI] -= 2 * h;
                                info.target->UpdateQuadPointCache();
                                Float leftSsVal = BdVal();

                                // Right, f(... x_i + h ... x_j - h ...)
                                perturbVertexI[indexI] += 2 * h;
                                perturbVertexJ[indexJ] -= 2 * h;
                                info.target->UpdateQuadPointCache();
                                Float rightSsVal = BdVal();

                                // Down, f(... x_i - h ... x_j - h ...)
                                perturbVertexI[indexI] -= 2 * h;
                                info.target->UpdateQuadPointCache();
                                Float downSsVal = BdVal();

                                Float hessEntry = (upSsVal - rightSsVal - leftSsVal + downSsVal) / (4 * hSq);
                                hess(matIndexI, matIndexJ) = hessEntry;
                            }

            return hess;
        }
#endif

    }
}