#include "bspch.h"
#include "BSSurface.h"

#include "Solver/Solver.h"
#include "Solver/SolverConfig.h"
#include "Utility/Utility.h"

namespace BSIPC
{
    BSSurface::BSSurface(SolverConfig config, UInt index, const BSSurfaceInfo& info)
    {
        std::vector<BSIPC::Vec2> control;               // material space
        std::vector<BSIPC::Vec3> initControl;           // world space

        Float clothWidth = info.restSizeX;
        Float clothHeight = info.restSizeY;
        UInt clothResX = info.resolutionX;
        UInt clothResY = info.resolutionY;

        if (info.irregularDomain)
        {
            for (const Vec2 vec : info.restNodePos)
            {
                control.push_back(vec);

                Vec3 control = Vec3(info.initRatioX * vec[0], info.initRatioY * vec[1], 0) + info.restPos;
                for (const auto& transformation : info.transformations)
                    transformation.ApplyTo(control);

                initControl.push_back(control);
            }
        }
        else
        {
            // Init material space control points
            for (UInt j = 0; j < clothResY; j++)
                for (UInt i = 0; i < clothResX; i++)
                {
                    Vec2 pos = Vec2(i * clothWidth / (clothResX - 1), j * clothHeight / (clothResY - 1));
                    control.push_back(pos);
                }

            for (auto& vec : control)
            {
                Vec3 control = Vec3(info.initRatioX * vec[0], info.initRatioY * vec[1], 0) + info.restPos;
                for (const auto& transformation : info.transformations)
                    transformation.ApplyTo(control);

                initControl.push_back(control);
            }
        }

        if (!info.startNodePos.empty())
        {
            BSIPC_ASSERT(initControl.size() == info.startNodePos.size(), "specified node position size and #{node number} mismatch");
            initControl = info.startNodePos;

            for (Vec3& pos : initControl)
                for (const auto& transformation : info.transformations)
                    transformation.ApplyTo(pos);
        }

        String sizeStr = fmt::format("material space/world space/surface info/ have inconsistent size: {}/{}/{}", 
            control.size(), initControl.size(), info.resolutionX * info.resolutionY);

        BSIPC_ASSERT(control.size() == initControl.size() && initControl.size() == info.resolutionX * info.resolutionY, sizeStr);

        this->irregular = info.irregularDomain;
        this->resolutionX = clothResX;
        this->resolutionY = clothResY;
        this->linearResolutionX = info.linearResolutionX;
        this->linearResolutionY = info.linearResolutionY;
        this->surfaceArea = info.surfaceArea;                     // Later this should be computed upon initializing BSSurface
        this->restControlPoints = control;
        this->controlPoints = initControl;
        this->solver = nullptr;
    }

    Vec3 BSSurface::EvalAt(Float u, Float v) const
    {
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        // Domain Start Points
        Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionX - 3));
        Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionY - 3));
        std::array<UInt, 3> uOffset = this->SelectSupportedOffset(u);
        std::array<UInt, 3> vOffset = this->SelectSupportedOffset(v);
        
        // Evaluate B-Spline 
        std::array<BSBasisFn, 3> horz = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasis::BasisFns);
        std::array<BSBasisFn, 3> vert = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasis::BasisFns);

        Vec3 result = Vec3(0., 0., 0.);

        for (Int i = 0; i < 3; i++)
        {
            for (Int j = 0; j < 3; j++)
            {
                Float uEval = (u - uFloor) + uOffset[i];
                Float vEval = (v - vFloor) + vOffset[j];

                Vec3 controlPoint = this->GetControlPoint(uFloor + i, vFloor + j);

                Float horzVal = horz[i](uEval);
                Float vertVal = vert[j](vEval);
                result += controlPoint * horzVal * vertVal;
            }
        }

        return result;
    }

    Vec2 BSSurface::EvalAtRest(Float u, Float v) const
    {
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        // Domain Start Points
        Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionX - 3));
        Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionY - 3));
        std::array<UInt, 3> uOffset = this->SelectSupportedOffset(u);
        std::array<UInt, 3> vOffset = this->SelectSupportedOffset(v);

        // Evaluate B-Spline 
        std::array<BSBasisFn, 3> horz = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasis::BasisFns);
        std::array<BSBasisFn, 3> vert = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasis::BasisFns);

        Vec2 result = Vec2(0., 0.);

        for (Int i = 0; i < 3; i++)
        {
            for (Int j = 0; j < 3; j++)
            {
                Float uEval = (u - uFloor) + uOffset[i];
                Float vEval = (v - vFloor) + vOffset[j];

                Vec2 controlPoint = this->GetRestControlPoint(uFloor + i, vFloor + j);

                Float horzVal = horz[i](uEval);
                Float vertVal = vert[j](vEval);

                result += controlPoint * horzVal * vertVal;
            }
        }

        return result;
    }

    Float BSSurface::EvalCoeffAt(UInt i, UInt j, Float u, Float v) const
    {
        if (!CoordInNodeClosedSupport(i, j, u, v))
            return 0.;

        // Domain Start Points
        Vec2 uSupport = this->SingleDimSupport(i, this->resolutionX);
        Vec2 vSupport = this->SingleDimSupport(j, this->resolutionY);

        // Evaluate B-Spline
        BSBasisFn horz = this->SelectBasisFnByIndex(i, BSBasis::BasisFns, this->resolutionX);
        BSBasisFn vert = this->SelectBasisFnByIndex(j, BSBasis::BasisFns, this->resolutionY);

        Float uEval = u - uSupport[0];
        Float vEval = v - vSupport[0];

        Float result = horz(uEval) * vert(vEval);

        return result;
    }

    Float BSSurface::EvalCoeffAt(UInt index, Float u, Float v) const
    {
        IVec2 ij = this->GetControlPointPosIndex(index);
        return this->EvalCoeffAt(ij[0], ij[1], u, v);
    }

    Vec3 BSSurface::NormalAt(Float u, Float v) const
    {
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        // Domain Start Points
        Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionX - 3));
        Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionY - 3));
        std::array<UInt, 3> uOffset = this->SelectSupportedOffset(u);
        std::array<UInt, 3> vOffset = this->SelectSupportedOffset(v);

        // Evaluate B-Spline 
        std::array<BSBasisFn, 3> horz = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasis::BasisFns);
        std::array<BSBasisFn, 3> vert = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasis::BasisFns);
        std::array<BSBasisFn, 3> horzDeriv = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasisDerivative::BasisFns);
        std::array<BSBasisFn, 3> vertDeriv = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasisDerivative::BasisFns);

        Vec3 uDeriv = Vec3(0., 0., 0.), vDeriv = Vec3(0., 0., 0.);

        for (Int i = 0; i < 3; i++)
        {
            for (Int j = 0; j < 3; j++)
            {
                Float uEval = (u - uFloor) + uOffset[i];
                Float vEval = (v - vFloor) + vOffset[j];

                Vec3 controlPoint = this->GetControlPoint(uFloor + i, vFloor + j);
                Float horzVal = horz[i](uEval);
                Float vertVal = vert[j](vEval);
                Float horzDerivVal = horzDeriv[i](uEval);
                Float vertDerivVal = vertDeriv[j](vEval);

                uDeriv += controlPoint * horzDerivVal * vertVal;
                vDeriv += controlPoint * horzVal * vertDerivVal;
            }
        }

        Vec3 normal = uDeriv.cross(vDeriv);
        normal.normalize();

        return normal;
    }

    Float BSSurface::EvalDistancebetween(Float u, Float v, Vec2 pos) const
    {
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        Vec2 testPos = EvalAtRest(u, v);
        return (testPos - pos).norm();
    }

    Vec2 BSSurface::PD_Distance_ParamCoord(Float u, Float v, Vec2 pos) const
    {
        if (!(u >= 0 && u <= this->UDomainMax()))
            BSIPC_INFO("u: {}, umax: {}", u, this->UDomainMax());

        if (!(v >= 0 && v <= this->VDomainMax()))
            BSIPC_INFO("v: {}, vmax: {}", v, this->VDomainMax());

        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        Vec2 testPos = EvalAtRest(u, v);

        // Domain Start Points
        Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionX - 3));
        Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionY - 3));
        std::array<UInt, 3> uOffset = this->SelectSupportedOffset(u);
        std::array<UInt, 3> vOffset = this->SelectSupportedOffset(v);

        // Evaluate B-Spline 
        std::array<BSBasisFn, 3> horz = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasis::BasisFns);
        std::array<BSBasisFn, 3> vert = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasis::BasisFns);
        std::array<BSBasisFn, 3> horzDeriv = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasisDerivative::BasisFns);
        std::array<BSBasisFn, 3> vertDeriv = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasisDerivative::BasisFns);

        // First Calculate coefficient term
        Vec2 PD_UCoeff = Vec2(0., 0.), PD_VCoeff = Vec2(0., 0.);

        for (Int i = 0; i < 3; i++)
        {
            for (Int j = 0; j < 3; j++)
            {
                Float uEval = (u - uFloor) + uOffset[i];
                Float vEval = (v - vFloor) + vOffset[j];

                Vec3 controlPoint = this->GetControlPoint(uFloor + i, vFloor + j);
                Float horzVal = horz[i](uEval);
                Float vertVal = vert[j](vEval);
                Float horzDerivVal = horzDeriv[i](uEval);
                Float vertDerivVal = vertDeriv[j](vEval);

                PD_UCoeff += this->GetRestControlPoint(uFloor + i, vFloor + j) * horzDerivVal * vertVal;
                PD_VCoeff += this->GetRestControlPoint(uFloor + i, vFloor + j) * horzVal * vertDerivVal;
            }
        }

        // Then compose
        Float drdu = 2 * (testPos[0] - pos[0]) * PD_UCoeff[0] + 2 * (testPos[1] - pos[1]) * PD_UCoeff[1];
        Float drdv = 2 * (testPos[0] - pos[0]) * PD_VCoeff[0] + 2 * (testPos[1] - pos[1]) * PD_VCoeff[1];

        return Vec2(drdu, drdv);
    }

    Vec2 BSSurface::ParametricCoordAtRest(Vec2 pos, Vec2 suggestUV) const
    {
        const Float tolerance = 1e-6;
        const UInt maxIter = 1000000;

        Float uCurrent = suggestUV[0];
        Float vCurrent = suggestUV[1];

        Float uTest = 0., vTest = 0.;
        Float alpha = 1.;

        for (UInt i = 0; i != maxIter; ++i)
        {
            Vec2 searchDir = this->PD_Distance_ParamCoord(uCurrent, vCurrent, pos);
            alpha = 2.;

            do
            {
                alpha *= 0.5;

                uTest = uCurrent - alpha * searchDir[0];
                vTest = vCurrent - alpha * searchDir[1];

                uTest = std::clamp(uTest, 0., static_cast<Float>(this->UDomainMax()));
                vTest = std::clamp(vTest, 0., static_cast<Float>(this->VDomainMax()));
            } while (this->EvalDistancebetween(uTest, vTest, pos) > this->EvalDistancebetween(uCurrent, vCurrent, pos));

            uCurrent = uTest;
            vCurrent = vTest;

            if (this->EvalDistancebetween(uTest, vTest, pos) < tolerance)
            {
                BSIPC_INFO("material to param coord #iter: {}", i);
                break;
            }

            if (i == maxIter - 1)
            {
                BSIPC_WARN("material to param coord reach max iter. Falling back to default (suggested) u-v coordinate");
                uCurrent = suggestUV[0];
                vCurrent = suggestUV[1];
            }
        }

        return Vec2(uCurrent, vCurrent);
    }

    Mat<3, 2> BSSurface::DeformationGradientAt(Float u, Float v, const Mat<2, 2>& pd_Param_Material)
    {
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        Float dudX1 = pd_Param_Material(0, 0), dudX2 = pd_Param_Material(0, 1), dvdX1 = pd_Param_Material(1, 0), dvdX2 = pd_Param_Material(1, 1);

        // TODO avoid repeating this in PDParamCoord and Eval
        Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionX - 3));
        Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionY - 3));
        auto uOffset = this->SelectSupportedOffset(u);
        auto vOffset = this->SelectSupportedOffset(v);

        // Basis Functions Used
        auto horzDeriv = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasisDerivative::BasisFns);
        auto vertDeriv = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasisDerivative::BasisFns);
        auto horzBasis = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasis::BasisFns);
        auto vertBasis = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasis::BasisFns);

        Float F11 = 0, F12 = 0, F21 = 0, F22 = 0, F31 = 0, F32 = 0;

        for (Int i = 0; i < 3; i++)
        {
            for (Int j = 0; j < 3; j++)
            {
                Float uEval = (u - uFloor) + uOffset[i];
                Float vEval = (v - vFloor) + vOffset[j];

                Vec3 controlPoint = this->GetControlPoint(uFloor + i, vFloor + j);
                Float horzVal = horzBasis[i](uEval);
                Float vertVal = vertBasis[j](vEval);
                Float horzDerivVal = horzDeriv[i](uEval);
                Float vertDerivVal = vertDeriv[j](vEval);

                F11 += horzDerivVal * vertVal * controlPoint[0] * dudX1 + horzVal * vertDerivVal * controlPoint[0] * dvdX1;
                F12 += horzDerivVal * vertVal * controlPoint[0] * dudX2 + horzVal * vertDerivVal * controlPoint[0] * dvdX2;
                F21 += horzDerivVal * vertVal * controlPoint[1] * dudX1 + horzVal * vertDerivVal * controlPoint[1] * dvdX1;
                F22 += horzDerivVal * vertVal * controlPoint[1] * dudX2 + horzVal * vertDerivVal * controlPoint[1] * dvdX2;
                F31 += horzDerivVal * vertVal * controlPoint[2] * dudX1 + horzVal * vertDerivVal * controlPoint[2] * dvdX1;
                F32 += horzDerivVal * vertVal * controlPoint[2] * dudX2 + horzVal * vertDerivVal * controlPoint[2] * dvdX2;
            }
        }

        Mat<3, 2> deformationGradient;
        deformationGradient <<
            F11, F12,
            F21, F22,
            F31, F32;

#if defined BSIPC_NUMERIC_TEST
        //Float du = .001;
        //if (u + du <= this->UDomainMax() && v + du <= this->UDomainMax())
        //{
        //    BSIPC_INFO("analytic: {}\n", ToStr(deformationGradient));
        //    Vec2 curAtRest = this->EvalAtRest(u, v);
        //    Vec2 uTestAtRest = this->EvalAtRest(u + du, v);
        //    Vec2 vTestAtRest = this->EvalAtRest(u, v + du);
        //    Float dx = uTestAtRest[0] - curAtRest[0];
        //    Float dy = vTestAtRest[1] - curAtRest[1];

        //    Vec3 cur = this->EvalAt(u, v);
        //    Vec3 uTest = this->EvalAt(u + du, v);
        //    Vec3 vTest = this->EvalAt(u, v + du);
        //    Vec3 dxVec = uTest - cur;
        //    Vec3 dyVec = vTest - cur;

        //    Mat<3, 2> deformationGradientTest;
        //    deformationGradientTest.col(0) = dxVec / dx;
        //    deformationGradientTest.col(1) = dyVec / dy;

        //    BSIPC_INFO("numerical: {}\n", ToStr(deformationGradientTest));
        //}
#endif

        return deformationGradient;
    }

    std::array<BSBasisFn, 3> BSSurface::SelectSupportedBasisFn(Float u, UInt bound, std::array<BSBasisFn, 5> basis) const
    {
        UInt uFloor = std::clamp(static_cast<UInt>(std::floor(u)), std::numeric_limits<UInt>::min(), static_cast<UInt>(bound - 3));
        std::array<BSBasisFn, 3> result;

        result[0] = basis[2];
        result[1] = basis[2];
        result[2] = basis[2];
        if (uFloor == 0)
        {
            result[0] = basis[0];
            result[1] = basis[1];
        }
        else if (uFloor == 1)
            result[0] = basis[1];
        else if (uFloor == bound - 3 || uFloor == bound - 2)
        {
            result[1] = basis[3];
            result[2] = basis[4];
        }
        else if (uFloor == bound - 4)
            result[2] = basis[3];

        return result;
    }

    std::array<UInt, 3> BSSurface::SelectSupportedOffset(Float u) const
    {
        std::array<UInt, 3> result = { 2, 1, 0 };
        UInt uFloor = static_cast<UInt>(std::floor(u));

        if (uFloor == 0)
            result = { 0, 0, 0 };
        else if (uFloor == 1)
            result = { 1, 1, 0 };

        return result;
    }

    Mat<2, 2> BSSurface::PD_Param_Material(Float u, Float v) const
    {
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound " + std::to_string(u) + " " + std::to_string(this->UDomainMax()));
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound " + std::to_string(v) + " " + std::to_string(this->VDomainMax()));

        Mat<2, 2> jacobian = this->PD_Material_Param(u, v);
        Mat<2, 2> invJacobian = jacobian.inverse();

        return invJacobian;
    }

    Mat<2, 3> BSSurface::PD2_Param_Material(
        Mat<2, 2> pd_Param_Material,
        Mat<2, 2> pd_Material_Param,
        Mat<2, 3> pd2_Material_Param
    ) const
    {
        // Known partial derivatives
        Float ux = pd_Param_Material(0, 0), uy = pd_Param_Material(0, 1), vx = pd_Param_Material(1, 0), vy = pd_Param_Material(1, 1);
        Float xu = pd_Material_Param(0, 0), xv = pd_Material_Param(1, 0), yu = pd_Material_Param(0, 1), yv = pd_Material_Param(1, 1);
        Float xuu = pd2_Material_Param(0, 0), xuv = pd2_Material_Param(0, 1), xvv = pd2_Material_Param(0, 2),
            yuu = pd2_Material_Param(1, 0), yuv = pd2_Material_Param(1, 1), yvv = pd2_Material_Param(1, 2);

        Float divisor = (xv * yu - xu * yv) * (xv * yu - xu * yv);

        Float xxxCoeff = xuu * yv * yv - 2 * xuv * yv * yu + xvv * yu * yu;
        Float xxyCoeff = yuu * yv * yv - 2 * yuv * yv * yu + yvv * yu * yu;
        Float xyxCoeff = xu * (xuv * yv - xvv * yu) + xv * (xuv * yu - xuu * yv);
        Float xyyCoeff = xu * (yuv * yv - yvv * yu) + xv * (yuv * yu - yuu * yv);
        Float yyxCoeff = xuu * xv * xv - 2 * xuv * xv * xu + xvv * xu * xu;
        Float yyyCoeff = yvv * xu * xu - 2 * yuv * xu * xv + yuu * xv * xv;

        Float uxx = -(ux * xxxCoeff + uy * xxyCoeff) / divisor;
        Float uxy = -(ux * xyxCoeff + uy * xyyCoeff) / divisor;
        Float uyy = -(ux * yyxCoeff + uy * yyyCoeff) / divisor;
        Float vxx = -(vx * xxxCoeff + vy * xxyCoeff) / divisor;
        Float vxy = -(vx * xyxCoeff + vy * xyyCoeff) / divisor;
        Float vyy = -(vx * yyxCoeff + vy * yyyCoeff) / divisor;

        Mat<2, 3> result;
        result << 
            uxx, uxy, uyy,
            vxx, vxy, vyy;
        return result;
    }

    Mat<2, 2> BSSurface::PD_Material_Param(Float u, Float v) const
    {
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        // Domain Start Points
        Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionX - 3));
        Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionY - 3));
        auto uOffset = this->SelectSupportedOffset(u);
        auto vOffset = this->SelectSupportedOffset(v);

        // Basis Functions Used
        auto horzDerivatives = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasisDerivative::BasisFns);
        auto vertDerivatives = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasisDerivative::BasisFns);
        auto horzBasis = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasis::BasisFns);
        auto vertBasis = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasis::BasisFns);

        Float dX1du = 0, dX1dv = 0, dX2du = 0, dX2dv = 0;

        for (Int i = 0; i < 3; i++)
        {
            for (Int j = 0; j < 3; j++)
            {
                Float uEval = (u - uFloor) + uOffset[i];
                Float vEval = (v - vFloor) + vOffset[j];

                Vec2 restControlPoint = this->GetRestControlPoint(uFloor + i, vFloor + j);
                Float horzVal = horzBasis[i](uEval);
                Float vertVal = vertBasis[j](vEval);
                Float horzDeriv = horzDerivatives[i](uEval);
                Float vertDeriv = vertDerivatives[j](vEval);

                dX1du += horzDeriv * vertVal * restControlPoint[0];
                dX1dv += horzVal * vertDeriv * restControlPoint[0];
                dX2du += horzDeriv * vertVal * restControlPoint[1];
                dX2dv += horzVal * vertDeriv * restControlPoint[1];
            }
        }

        Mat<2, 2> jacobian;
        jacobian <<
            dX1du, dX1dv,
            dX2du, dX2dv;

        return jacobian;
    }

    Mat<2, 3> BSSurface::PD2_Material_Param(Float u, Float v) const
    {
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        // Domain Start Points
        Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionX - 3));
        Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionY - 3));
        auto uOffset = this->SelectSupportedOffset(u);
        auto vOffset = this->SelectSupportedOffset(v);

        // Basis Functions Used
        auto horzBasisFn = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasis::BasisFns);
        auto vertBasisFn = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasis::BasisFns);
        auto horzDerivFn = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasisDerivative::BasisFns);
        auto vertDerivFn = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasisDerivative::BasisFns);
        auto horzSecondDerivFn = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasisSecondDerivative::BasisFns);
        auto vertSecondDerivFn = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasisSecondDerivative::BasisFns);

        Float d2X1du2 = 0, d2X1dv2 = 0, d2X2du2 = 0, d2X2dv2 = 0, d2X1dudv = 0, d2X2dudv = 0;

        for (UInt i = 0; i != 3; ++i)
        {
            for (UInt j = 0; j != 3; ++j)
            {
                Float uEval = (u - uFloor) + uOffset[i];
                Float vEval = (v - vFloor) + vOffset[j];

                Vec2 restControlPoint = this->GetRestControlPoint(uFloor + i, vFloor + j);
                Float horzVal = horzBasisFn[i](uEval);
                Float vertVal = vertBasisFn[j](vEval);
                Float horzDeriv = horzDerivFn[i](uEval);
                Float vertDeriv = vertDerivFn[j](vEval);
                Float horzSecondDeriv = horzSecondDerivFn[i](uEval);
                Float vertSecondDeriv = vertSecondDerivFn[j](vEval);

                d2X1du2 += horzSecondDeriv * vertVal * restControlPoint[0];
                d2X1dv2 += horzVal * vertSecondDeriv * restControlPoint[0];
                d2X1dudv += horzDeriv * vertDeriv * restControlPoint[0];
                d2X2dudv += horzDeriv * vertDeriv * restControlPoint[1];
                d2X2du2 += horzSecondDeriv * vertVal * restControlPoint[1];
                d2X2dv2 += horzVal * vertSecondDeriv * restControlPoint[1];
            }
        }

        Mat<2, 3> result;
        result <<
            d2X1du2, d2X1dudv, d2X1dv2,
            d2X2du2, d2X2dudv, d2X2dv2;

        return result;
    }

    Mat<3, 2> BSSurface::PD_World_Param(Float u, Float v) const
    {
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        // Domain Start Points
        Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionX - 3));
        Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionY - 3));
        auto uOffset = this->SelectSupportedOffset(u);
        auto vOffset = this->SelectSupportedOffset(v);

        // Basis Functions Used
        auto horzDerivatives = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasisDerivative::BasisFns);
        auto vertDerivatives = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasisDerivative::BasisFns);
        auto horzBasis = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasis::BasisFns);
        auto vertBasis = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasis::BasisFns);

        Float dxdu = 0, dxdv = 0, dydu = 0, dydv = 0, dzdu = 0, dzdv = 0;

        for (UInt i = 0; i != 3; ++i)
        {
            for (UInt j = 0; j != 3; ++j)
            {
                Float uEval = (u - uFloor) + uOffset[i];
                Float vEval = (v - vFloor) + vOffset[j];

                Vec3 controlPoint = this->GetControlPoint(uFloor + i, vFloor + j);
                Float horzVal = horzBasis[i](uEval);
                Float vertVal = vertBasis[j](vEval);
                Float horzDeriv = horzDerivatives[i](uEval);
                Float vertDeriv = vertDerivatives[j](vEval);

                dxdu += horzDeriv * vertVal * controlPoint[0];
                dxdv += horzVal * vertDeriv * controlPoint[0];
                dydu += horzDeriv * vertVal * controlPoint[1];
                dydv += horzVal * vertDeriv * controlPoint[1];
                dzdu += horzDeriv * vertVal * controlPoint[2];
                dzdv += horzVal * vertDeriv * controlPoint[2];
            }
        }

        Mat<3, 2> jacobian;
        jacobian <<
            dxdu, dxdv,
            dydu, dydv,
            dzdu, dzdv;

        return jacobian;
    }

    Mat<3, 3> BSSurface::PD2_World_Param(Float u, Float v) const
    {
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        // Domain Start Points
        Int uFloor = std::clamp(static_cast<Int>(std::floor(u)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionX - 3));
        Int vFloor = std::clamp(static_cast<Int>(std::floor(v)), std::numeric_limits<Int>::min(), static_cast<Int>(this->resolutionY - 3));
        auto uOffset = this->SelectSupportedOffset(u);
        auto vOffset = this->SelectSupportedOffset(v);

        // Basis Functions Used
        auto horzBasisFn = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasis::BasisFns);
        auto vertBasisFn = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasis::BasisFns);
        auto horzDerivFn = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasisDerivative::BasisFns);
        auto vertDerivFn = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasisDerivative::BasisFns);
        auto horzSecondDerivFn = this->SelectSupportedBasisFn(u, this->resolutionX, BSBasisSecondDerivative::BasisFns);
        auto vertSecondDerivFn = this->SelectSupportedBasisFn(v, this->resolutionY, BSBasisSecondDerivative::BasisFns);

        Float d2xdu2 = 0, d2xdv2 = 0, d2xdudv = 0, d2ydu2 = 0, d2ydv2 = 0, d2ydudv = 0, d2zdu2 = 0, d2zdv2 = 0, d2zdudv = 0;

        for (UInt i = 0; i != 3; ++i)
        {
            for (UInt j = 0; j != 3; ++j)
            {
                Float uEval = (u - uFloor) + uOffset[i];
                Float vEval = (v - vFloor) + vOffset[j];

                Vec3 controlPoint = this->GetControlPoint(uFloor + i, vFloor + j);
                Float horzVal = horzBasisFn[i](uEval);
                Float vertVal = vertBasisFn[j](vEval);
                Float horzDeriv = horzDerivFn[i](uEval);
                Float vertDeriv = vertDerivFn[j](vEval);
                Float horzSecondDeriv = horzSecondDerivFn[i](uEval);
                Float vertSecondDeriv = vertSecondDerivFn[j](vEval);

                d2xdu2 += horzSecondDeriv * vertVal * controlPoint[0];
                d2xdv2 += horzVal * vertSecondDeriv * controlPoint[0];
                d2xdudv += horzDeriv * vertDeriv * controlPoint[0];
                d2ydu2 += horzSecondDeriv * vertVal * controlPoint[1];
                d2ydv2 += horzVal * vertSecondDeriv * controlPoint[1];
                d2ydudv += horzDeriv * vertDeriv * controlPoint[1];
                d2zdu2 += horzSecondDeriv * vertVal * controlPoint[2];
                d2zdv2 += horzVal * vertSecondDeriv * controlPoint[2];
                d2zdudv += horzDeriv * vertDeriv * controlPoint[2];
            }
        }

        Mat<3, 3> hessian;
        hessian <<
            d2xdu2, d2xdudv, d2xdv2,
            d2ydu2, d2ydudv, d2ydv2,
            d2zdu2, d2zdudv, d2zdv2;

        return hessian;
    }

    Vec2 BSSurface::PD_DeformGrad_CtrlPt_At(UInt i, UInt j, Float u, Float v) const
    {
        BSIPC_ASSERT(i >= 0 && i <= this->resolutionX - 1, "Control Point i out of bound");
        BSIPC_ASSERT(j >= 0 && j <= this->resolutionY - 1, "Control Point j out of bound");
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        if (!CoordInNodeClosedSupport(i, j, u, v))
            return Vec2(0., 0.);

        Mat<2, 2> jacobian = this->PD_Param_Material(u, v);
        Float dudX1 = jacobian(0, 0), dudX2 = jacobian(0, 1), dvdX1 = jacobian(1, 0), dvdX2 = jacobian(1, 1);

        // Domain Start Points
        Vec2 uSupport = SingleDimSupport(i, this->resolutionX);
        Vec2 vSupport = SingleDimSupport(j, this->resolutionY);

        Float uEval = u - uSupport[0];
        Float vEval = v - vSupport[0];

        // Select Basis Functions
        BSBasisFn horzFn = SelectBasisFnByIndex(i, BSBasis::BasisFns, this->resolutionX);
        BSBasisFn vertFn = SelectBasisFnByIndex(j, BSBasis::BasisFns, this->resolutionY);
        BSBasisFn horzDeriv = SelectBasisFnByIndex(i, BSBasisDerivative::BasisFns, this->resolutionX);
        BSBasisFn vertDeriv = SelectBasisFnByIndex(j, BSBasisDerivative::BasisFns, this->resolutionY);

        Float uCoeff = horzDeriv(uEval) * vertFn(vEval);
        Float vCoeff = horzFn(uEval) * vertDeriv(vEval);

        Float result1 = uCoeff * dudX1 + vCoeff * dvdX1;
        Float result2 = uCoeff * dudX2 + vCoeff * dvdX2;

        return Vec2(result1, result2);
    }

    void BSSurface::PrepareQuadPoints()
    {
        const SolverConfig& config = this->solver->GetConfig();
        UInt order = config.quadOrder;
        QuadScheme scheme = config.quadScheme;
        QuadNodesData quadData1D = GaussianQuadrature::GetNodesOfOrder(order);
        Float w = 0., patchW = 0.;
        UInt numQuadPoints = 0;

        tbb::mutex nearbyCacheMutex;

        UInt ctrlPtCnt = this->GetControlPointCnt();
        this->nearbyQuadPointCache.resize(ctrlPtCnt);

        // rectangular surface for quadrature on irregular domain
        BSSurfaceInfo quadSurfInfo;
        quadSurfInfo.resolutionX = this->resolutionX;
        quadSurfInfo.resolutionY = this->resolutionY;
        //quadSurfInfo.restSizeX = std::abs(this->EvalAtRest(uDomain, 0)[0] - this->EvalAtRest(0, 0)[0]);
        //quadSurfInfo.restSizeY = std::abs(this->EvalAtRest(0, vDomain)[1] - this->EvalAtRest(0, 0)[1]);
        quadSurfInfo.restSizeX = 1.;
        quadSurfInfo.restSizeY = 1.;

        BSSurface quadSurf(this->solver->GetConfig(), std::numeric_limits<UInt>::max(), quadSurfInfo);

        // Init higher resolution quad points
        if (IsQuadSchemeGlobal(scheme)) // Global
        {
            BSIPC_ASSERT(!this->irregular, "global quadrature scheme is not supported for irregular B-spline surfaces");

            Float domainWidth = this->LengthSideX();
            Float domainHeight = this->LengthSideY();

            BSIPC_INFO("Domain Width: {:02.2f}, Domain Height: {:02.2f}", domainWidth, domainHeight);

            // widthData and heightData stores quadrature points in parametric coordinates
            QuadNodesData widthData = quadData1D, heightData = quadData1D;

            // Transforming quad points from [-1, 1] to surface domain (in world coordinates)
            for (UInt i = 0; i != order; ++i)
            {
                widthData.nodes[i] = TransformSamplePosition(quadData1D.nodes[i], 0, domainWidth);
                //BSIPC_INFO("U Node (World): {:02.4f}", widthData.nodes[i]);
                Vec2 materialCoord = Vec2(widthData.nodes[i], 0);
                Vec2 suggestUV = Vec2(this->UDomainMax() * materialCoord[0] / this->LengthSideX(),
                    this->VDomainMax() * materialCoord[1] / this->LengthSideY());
                Vec2 paramCoord = this->ParametricCoordAtRest(materialCoord, suggestUV);
                widthData.nodes[i] = paramCoord[0];
            }

            for (UInt i = 0; i != order; ++i)
            {
                heightData.nodes[i] = TransformSamplePosition(quadData1D.nodes[i], 0, domainHeight);
                //BSIPC_INFO("V Node (World): {:02.4f}", heightData.nodes[i]);
                Vec2 materialCoord = Vec2(0, heightData.nodes[i]);
                Vec2 suggestUV = Vec2(this->UDomainMax() * materialCoord[0] / this->LengthSideX(),
                    this->VDomainMax() * materialCoord[1] / this->LengthSideY());
                Vec2 paramCoord = this->ParametricCoordAtRest(materialCoord, suggestUV);
                heightData.nodes[i] = paramCoord[1];
            }

            // Trim (Merge) quadrature points according to the config
            const UInt cellLimit = 2;
            const UInt cornerLimit = 3;
            if (scheme == QuadScheme::GLOBAL_UNIFORM)
                ;   
            else if (scheme == QuadScheme::GLOBAL_LIMITED)
            {
                // Do not need to fill in specLimits here as the limit is uniform
                std::map<UInt, UInt> specLimits;
                widthData = TrimDenseQuadPoints(specLimits, widthData, cellLimit, this->UDomainMax());
                heightData = TrimDenseQuadPoints(specLimits, heightData, cellLimit, this->VDomainMax());
            }
            else if (scheme == QuadScheme::GLOBAL_EDGE_DENSE)
            {
                std::map<UInt, UInt> specLimits;
                specLimits[0] = cornerLimit;
                specLimits[this->UDomainMax() - 1] = cornerLimit;
                widthData = TrimDenseQuadPoints(specLimits, widthData, cellLimit, this->UDomainMax());

                specLimits.clear();

                specLimits[0] = cornerLimit;
                specLimits[this->VDomainMax() - 1] = cornerLimit;
                heightData = TrimDenseQuadPoints(specLimits, heightData, cellLimit, this->VDomainMax());
            }
            else if (scheme == QuadScheme::GLOBAL_CHESSBOARD)
            {
                std::map<UInt, UInt> specLimits;
                for (UInt i = 0; i != this->UDomainMax(); ++i)
                    specLimits[i] = i % 2 == 0 ? 2 : 1;
                specLimits[0] = cornerLimit;
                specLimits[this->UDomainMax() - 1] = cornerLimit;
                widthData = TrimDenseQuadPoints(specLimits, widthData, cellLimit, this->UDomainMax());

                specLimits.clear();

                for (UInt i = 0; i != this->VDomainMax(); ++i)
                    specLimits[i] = i % 2 == 0 ? 1 : 2;
                specLimits[0] = cornerLimit;
                specLimits[this->VDomainMax() - 1] = cornerLimit;
                heightData = TrimDenseQuadPoints(specLimits, heightData, cellLimit, this->VDomainMax());
            }

            Float patchArea = domainWidth * domainHeight;
            Float patchWeightCoeff = patchArea / 4.0;

            // Convert quadrature points from world coordaintes to parametric coordinates for evaluation
            for (UInt i = 0; i != widthData.nodes.size(); ++i)
            {
                for (UInt j = 0; j != heightData.nodes.size(); ++j)
                {
                    Float curQuadWeight = widthData.weights[i] * heightData.weights[j] * patchWeightCoeff;
                    //BSIPC_INFO("(Parametric Coordinates) U: {:02.2f}, V: {:02.2f}, Weight: {:02.4f}",
                    //    widthData.nodes[i], heightData.nodes[j], curQuadWeight);
                    Mat<2, 2> pd_Material_Param = this->PD_Material_Param(widthData.nodes[i], heightData.nodes[j]);
                    Mat<2, 3> pd2_Material_Param = this->PD2_Material_Param(widthData.nodes[i], heightData.nodes[j]);
                    Mat<2, 2> pd_Param_Material = this->PD_Param_Material(widthData.nodes[i], heightData.nodes[j]);
                    Mat<2, 3> pd2_Param_Material = this->PD2_Param_Material(
                        pd_Param_Material, pd_Material_Param, pd2_Material_Param
                    );
                    this->ssQuadPoints.emplace_back(QuadPoint
                        {
                            widthData.nodes[i], heightData.nodes[j], curQuadWeight,
                            pd_Param_Material, pd2_Param_Material, pd_Material_Param, pd2_Material_Param
                        });

                    // push the quadrature point into the nearbyQuadCache to accelerate mass matrix initialization
                    UInt curQuadIndex = static_cast<UInt>(i * heightData.nodes.size() + j);
                    IVec2 supportStart = this->SupportStartOf(widthData.nodes[i], heightData.nodes[j]);

                    for (UInt k = 0; k != 9; ++k)   // for each control point whose support contains the quadrature point
                    {
                        UInt uOffset = k % 3;
                        UInt vOffset = k / 3;

                        IVec2 curSupport = supportStart + IVec2(uOffset, vOffset);
                        UInt curSupportIndex = this->GetControlPointIndex(curSupport[0], curSupport[1]);

                        BSIPC_ASSERT(
                            curSupportIndex < this->nearbyQuadPointCache.size(),
                            "Nearby Quad Point Cache size too small!");

                        this->nearbyQuadPointCache[curSupportIndex].insert(curQuadIndex);
                    }

                    w += curQuadWeight;
                    ++numQuadPoints;
                }
            }
        }
        else if (scheme == QuadScheme::LOCAL_DUAL_CHESSBOARD)
        {
            UInt uMax = this->UDomainMax();
            UInt vMax = this->VDomainMax();
            UInt patchSize = uMax * vMax;

            this->patchQuadStartIndices.resize(patchSize);
            BSIPC_WARN("{}", this->patchQuadStartIndices.size());

            UInt dualGridSizeX = uMax + 1, dualGridSizeY = vMax + 1;
            UInt iterSize = dualGridSizeX * dualGridSizeY;

            auto PatchIndexFromCoord = [uMax, vMax](Float u, Float v) -> Int
                {
                    BSIPC_ASSERT(u <= uMax && v <= vMax, "input outside domain of surface");

                    UInt uFloor = static_cast<UInt>(std::floor(u));
                    UInt vFloor = static_cast<UInt>(std::floor(v));

                    if (uFloor == uMax)
                        --uFloor;
                    if (vFloor == vMax)
                        --vFloor;

                    return vFloor * uMax + uFloor;
                };

            auto UVFromDualGridIndex = [dualGridSizeX, dualGridSizeY](UInt gridIndexI, UInt gridIndexJ) -> Vec2
                {
                    Vec2 result = Vec2::Zero();

                    if (gridIndexI == 1)
                        result[0] = 1;
                    else if (gridIndexI == dualGridSizeX - 1)
                        result[0] = gridIndexI - 1;
                    else
                        result[0] = gridIndexI - 0.5;

                    if (gridIndexJ == 1)
                        result[1] = 1;
                    else if (gridIndexJ == dualGridSizeY - 1)
                        result[1] = gridIndexJ - 1;
                    else
                        result[1] = gridIndexJ - 0.5;

                    return result;
                };

            // Init ssQuadPoints
            {
                // map: index of vector -> quadrature points on the surface
                // Vec3: [quadPt.u, quadPt.v, weight]
                std::vector<tbb::concurrent_vector<Vec3>> patchToQuadPoints(patchSize);

                Float w = 0;

                // Init patches on the boundary
                {
                    tbb::parallel_for(
                        0, static_cast<Int>(uMax * vMax), 1, [&](Int iter)
                        //for (Int iter = 0; iter != static_cast<Int>(iterSize); ++ iter)
                        {
                            UInt i = static_cast<UInt>(iter % uMax);
                            UInt j = static_cast<UInt>(iter / uMax);

                            // Set up Quadrature computing domain
                            Float uMin = static_cast<Float>(i);
                            Float uMax = static_cast<Float>(i + 1);
                            Float vMin = static_cast<Float>(j);
                            Float vMax = static_cast<Float>(j + 1);

                            Vec2 llCorner = this->EvalAtRest(uMin, vMin);
                            Vec2 luCorner = this->EvalAtRest(uMin, vMax);
                            Vec2 ulCorner = this->EvalAtRest(uMax, vMin);
                            Vec2 uuCorner = this->EvalAtRest(uMax, vMax);

                            Float patchArea = 0;
                            patchArea += trigArea(llCorner, luCorner, ulCorner);
                            patchArea += trigArea(luCorner, ulCorner, uuCorner);

                            Float patchWeightCoeff = patchArea / 4.0;

                            auto [orderI, orderJ] = this->LocalPatchOrder(QuadScheme::LOCAL_EDGE_DENSE, i, j);

                            if (orderI == 3 || orderJ == 3) // the patch is on the boundary
                            {
                                QuadNodesData quadData1DI = GaussianQuadrature::GetNodesOfOrder(orderI);
                                QuadNodesData quadData1DJ = GaussianQuadrature::GetNodesOfOrder(orderJ);

                                // Compute quadrature points and convert to parametric coordinates
                                for (UInt k = 0; k != orderI; ++k)
                                    for (UInt l = 0; l != orderJ; ++l)
                                    {
                                        Vec2 curQuadPointParamCoord;
                                        Float curQuadWeight;

                                        Float curQuadPointX = TransformSamplePosition(quadData1DI.nodes[k], llCorner[0], uuCorner[0]);
                                        Float curQuadPointY = TransformSamplePosition(quadData1DJ.nodes[l], llCorner[1], uuCorner[1]);
                                        curQuadWeight = quadData1DI.weights[k] * quadData1DJ.weights[l] * patchWeightCoeff;

                                        Vec2 curQuadPoint = Vec2(curQuadPointX, curQuadPointY);
                                        Vec2 suggestUV = Vec2(uMin + 0.5, vMin + 0.5);
                                        curQuadPointParamCoord = this->ParametricCoordAtRest(curQuadPoint, suggestUV);

                                        BSIPC_INFO("Quad Point: x: {:02.2f} y: {:02.2f} u: {:02.2f} v: {:02.2f} w: {:04.8f}",
                                            curQuadPointX, curQuadPointY, curQuadPointParamCoord[0], curQuadPointParamCoord[1], curQuadWeight);

                                        UInt patchIndex = PatchIndexFromCoord(curQuadPointParamCoord[0], curQuadPointParamCoord[1]);
                                        patchToQuadPoints[patchIndex].push_back(Vec3(curQuadPointParamCoord[0], curQuadPointParamCoord[1], curQuadWeight));
                                    }
                            }
                        }
                    );
                }

                QuadNodesData order1 = GaussianQuadrature::GetNodesOfOrder(1);
                QuadNodesData order2 = GaussianQuadrature::GetNodesOfOrder(2);
                QuadNodesData order3 = GaussianQuadrature::GetNodesOfOrder(3);
                QuadNodesData order4 = GaussianQuadrature::GetNodesOfOrder(4);

                // Iterate through all interior quad points, fill into corresponding patch index.
                tbb::parallel_for(
                    0, static_cast<Int>(iterSize), 1, [&](Int idx)
                    {
                        Int gridIndexI = idx % dualGridSizeX;
                        Int gridIndexJ = idx / dualGridSizeX;

                        if (gridIndexI == 0 || gridIndexI == dualGridSizeX - 1 || gridIndexJ == 0 || gridIndexJ == dualGridSizeY - 1)
                            return;

                        Vec2 uvMin = UVFromDualGridIndex(gridIndexI, gridIndexJ);
                        Vec2 uvMax = UVFromDualGridIndex(gridIndexI + 1, gridIndexJ + 1);

                        Vec2 llCorner = this->EvalAtRest(uvMin[0], uvMin[1]);
                        Vec2 luCorner = this->EvalAtRest(uvMin[0], uvMax[1]);
                        Vec2 ulCorner = this->EvalAtRest(uvMax[0], uvMin[1]);
                        Vec2 uuCorner = this->EvalAtRest(uvMax[0], uvMax[1]);

                        Float patchArea = 0;
                        patchArea += trigArea(llCorner, luCorner, ulCorner);
                        patchArea += trigArea(luCorner, ulCorner, uuCorner);

                        Float patchWeightCoeff = patchArea / 4.0;

                        if ((gridIndexI == 1 || gridIndexI == dualGridSizeX - 2) && (gridIndexJ == 1 || gridIndexJ == dualGridSizeY - 2))
                        {
                            Vec2 quadPos = (llCorner + uuCorner + luCorner + ulCorner) / 4.0;
                            Vec2 suggestUV = (uvMin + uvMax) / 2.0;

                            Vec2 quadUVPos = this->ParametricCoordAtRest(quadPos, suggestUV);
                            UInt patchIndex = PatchIndexFromCoord(quadUVPos[0], quadUVPos[1]);
                            Float weight = order1.weights[0] * order1.weights[0] * patchWeightCoeff;
                            patchToQuadPoints[patchIndex].push_back(Vec3(quadUVPos[0], quadUVPos[1], weight));
                        }
                        else if (gridIndexI == 1 || gridIndexI == dualGridSizeX - 2)
                        {
                            //Float quadU = (uvMin[0] + uvMax[0]) / 2.;
                            Float quadU = gridIndexI == 1 ? uvMin[0] : uvMax[0];
                            for (UInt idx = 0; idx != 2; ++idx)
                            {
                                Float quadV = TransformSamplePosition(order2.nodes[idx], 0, 1) + uvMin[1];
                                Float weight = order1.weights[0] * order2.weights[idx] * patchWeightCoeff;
                                Int patchIndex = PatchIndexFromCoord(quadU, quadV);
                                patchToQuadPoints[patchIndex].push_back(Vec3(quadU, quadV, weight));
                            }
                        }
                        else if (gridIndexJ == 1 || gridIndexJ == dualGridSizeY - 2)
                        {
                            //Float quadV = (uvMin[1] + uvMax[1]) / 2.0;
                            Float quadV = gridIndexJ == 1 ? uvMin[1] : uvMax[1];
                            for (UInt idx = 0; idx != 2; ++idx)
                            {
                                Float quadU = TransformSamplePosition(order2.nodes[idx], 0, 1) + uvMin[0];
                                Float weight = order2.weights[idx] * order1.weights[0] * patchWeightCoeff;
                                Int patchIndex = PatchIndexFromCoord(quadU, quadV);
                                patchToQuadPoints[patchIndex].push_back(Vec3(quadU, quadV, weight));
                            }
                        }
                        else
                        {
                            UInt indicator = gridIndexI + gridIndexJ;
                            if (indicator % 2 == 0) // 3 * 1
                            {
                                for (UInt idx = 0; idx != 2; ++idx)
                                {
                                    Float quadU = TransformSamplePosition(order2.nodes[idx], 0, 1) + uvMin[0];
                                    Float quadV = (uvMin[1] + uvMax[1]) / 2.0;
                                    Float weight = order2.weights[idx] * order1.weights[0] * patchWeightCoeff;
                                    Int patchIndex = PatchIndexFromCoord(quadU, quadV);
                                    patchToQuadPoints[patchIndex].push_back(Vec3(quadU, quadV, weight));
                                }
                            }
                            else // 1 * 3
                            {
                                for (UInt idx = 0; idx != 2; ++idx)
                                {
                                    Float quadU = (uvMin[0] + uvMax[0]) / 2.0;
                                    Float quadV = TransformSamplePosition(order2.nodes[idx], 0, 1) + uvMin[1];
                                    Float weight = order1.weights[0] * order2.weights[idx] * patchWeightCoeff;
                                    Int patchIndex = PatchIndexFromCoord(quadU, quadV);
                                    patchToQuadPoints[patchIndex].push_back(Vec3(quadU, quadV, weight));
                                }
                            }
                        }
                    }
                );

                // Record start index of all patches
                UInt startIndex = 0;
                for (UInt i = 0; i != patchSize; ++i)
                {
                    this->patchQuadStartIndices[i] = startIndex;
                    startIndex += static_cast<UInt>(patchToQuadPoints[i].size());
                }
                // here startIndex = #quadPoints
                this->ssQuadPoints.resize(startIndex, QuadPoint::Zero());

                // Concurrently fill into the final QuadPoint vector
                //tbb::parallel_for(
                    //0, static_cast<Int>(patchToQuadPoints.size()), 1, [&](UInt idx)
                    for (UInt idx = 0; idx != patchToQuadPoints.size(); ++idx)
                    {
                        for (const Vec3& vec : patchToQuadPoints[idx])
                            BSIPC_INFO("u: {:02.2f} v: {:02.2f} w: {:04.8f}", vec[0], vec[1], vec[2]);

                        BSIPC_INFO_SEP

                        for (UInt j = 0; j != patchToQuadPoints[idx].size(); ++j)
                        {
                            Vec2 quadUVPos = Vec2(patchToQuadPoints[idx][j][0], patchToQuadPoints[idx][j][1]);
                            Float weight = patchToQuadPoints[idx][j][2];

                            Mat<2, 2> pd_Material_Param = this->PD_Material_Param(quadUVPos[0], quadUVPos[1]);
                            Mat<2, 3> pd2_Material_Param = this->PD2_Material_Param(quadUVPos[0], quadUVPos[1]);
                            Mat<2, 2> pd_Param_Material = this->PD_Param_Material(quadUVPos[0], quadUVPos[1]);
                            Mat<2, 3> pd2_Param_Material = this->PD2_Param_Material(
                                pd_Param_Material, pd_Material_Param, pd2_Material_Param
                            );

                            BSIPC_ASSERT(this->ssQuadPoints.size() < std::numeric_limits<UInt>::max(),
                                "Quadrature points size too large!");

                            UInt curQuadIndex = this->patchQuadStartIndices[idx] + j;
                            QuadPoint curQuad = QuadPoint
                                {
                                    quadUVPos[0], quadUVPos[1], weight,
                                    pd_Param_Material, pd2_Param_Material, pd_Material_Param, pd2_Material_Param
                                };

                            w += weight;

                            this->ssQuadPoints[curQuadIndex] = curQuad;
                            IVec2 supportStart = this->SupportStartOf(curQuad.U(), curQuad.V());
                            for (UInt n = 0; n != 9; ++n)   // for each control point whose support contains the quadrature point
                            {
                                UInt uOffset = n % 3;
                                UInt vOffset = n / 3;

                                IVec2 curSupport = supportStart + IVec2(uOffset, vOffset);
                                UInt curSupportIndex = this->GetControlPointIndex(curSupport[0], curSupport[1]);

                                BSIPC_ASSERT(
                                    curSupportIndex < this->nearbyQuadPointCache.size(),
                                    "Nearby Quad Point Cache size too small!");

                                this->nearbyQuadPointCache[curSupportIndex].insert(curQuadIndex);
                            }
                        }
                    }
                //);

                BSIPC_PEEK(w);
            }

            //for (const QuadPoint& quad : this->ssQuadPoints)
            //    BSIPC_INFO("u: {:02.2f} v: {:02.2f} w: {:04.8f}", quad.U(), quad.V(), quad.Weight());

            // Init bdQuadPoints
        }
        else // Local
        {
            // TODO fill local quad points in parallel
            BSIPC_INFO("Per-patch Quadrature Order: {}", order);

            UInt totalQuadPtCnt = this->QuadPointCnt(scheme);
            this->ssQuadPoints.resize(totalQuadPtCnt, QuadPoint::Zero());
            BSIPC_PEEK(totalQuadPtCnt);

            UInt uMax = this->UDomainMax();
            UInt vMax = this->VDomainMax();
            UInt iterSize = uMax * vMax;
            BSIPC_INFO("u: {}, v: {}", uMax, vMax);
            this->patchQuadStartIndices.resize(iterSize);

            tbb::parallel_for(
                0, static_cast<Int>(iterSize), 1, [&](Int iter) 
            //for (Int iter = 0; iter != static_cast<Int>(iterSize); ++ iter)
                {
                    UInt i = static_cast<UInt>(iter % uMax);
                    UInt j = static_cast<UInt>(iter / uMax);

                    UInt quadPtOffset = this->QuadPointStartIndexAtPatch(scheme, i, j);
                    this->patchQuadStartIndices[iter] = quadPtOffset;

                    // Set up Quadrature computing domain
                    Float uMin = static_cast<Float>(i);
                    Float uMax = static_cast<Float>(i + 1);
                    Float vMin = static_cast<Float>(j);
                    Float vMax = static_cast<Float>(j + 1);

                    Vec2 llCorner = this->EvalAtRest(uMin, vMin);
                    Vec2 luCorner = this->EvalAtRest(uMin, vMax);
                    Vec2 ulCorner = this->EvalAtRest(uMax, vMin);
                    Vec2 uuCorner = this->EvalAtRest(uMax, vMax);

                    //BSIPC_INFO("ll: {:02.4f}, {:02.4f}; lu: {:02.4f}, {:02.4f}; ul: {:02.4f}, {:02.4f}; uu: {:02.4f}, {:02.4f}", 
                    //    llCorner[0], llCorner[1],
                    //    luCorner[0], luCorner[1],
                    //    ulCorner[0], ulCorner[1],
                    //    uuCorner[0], uuCorner[1]
                    //    );

                    Float patchArea = 0;
                    patchArea += trigArea(llCorner, luCorner, ulCorner);
                    patchArea += trigArea(luCorner, ulCorner, uuCorner);

                    Float patchWeightCoeff = patchArea / 4.0;

                    auto [orderI, orderJ] = this->LocalPatchOrder(scheme, i, j);

                    QuadNodesData quadData1DI = GaussianQuadrature::GetNodesOfOrder(orderI);
                    QuadNodesData quadData1DJ = GaussianQuadrature::GetNodesOfOrder(orderJ);

                    // Compute quadrature points and convert to parametric coordinates
                    for (UInt k = 0; k != orderI; ++k)
                        for (UInt l = 0; l != orderJ; ++l)
                        {
                            Vec2 curQuadPointParamCoord;
                            Float curQuadWeight;

                            if (this->irregular)
                            {
                                Vec2 llQuadCorner = quadSurf.EvalAtRest(uMin, vMin);
                                Vec2 uuQuadCorner = quadSurf.EvalAtRest(uMax, vMax);

                                Float curQuadPointX = TransformSamplePosition(quadData1DI.nodes[k], llQuadCorner[0], uuQuadCorner[0]);
                                Float curQuadPointY = TransformSamplePosition(quadData1DJ.nodes[l], llQuadCorner[1], uuQuadCorner[1]);
                                curQuadWeight = quadData1DI.weights[k] * quadData1DJ.weights[l] * patchWeightCoeff;

                                Vec2 curQuadPoint = Vec2(curQuadPointX, curQuadPointY);
                                Vec2 suggestUV = Vec2(uMin + 0.5, vMin + 0.5);
                                BSIPC_INFO("from quad surface");
                                curQuadPointParamCoord = quadSurf.ParametricCoordAtRest(curQuadPoint, suggestUV);

                                BSIPC_INFO("Quad Point: x: {:02.2f} y: {:02.2f} u: {:02.2f} v: {:02.2f} w: {:04.8f}",
                                    curQuadPointX, curQuadPointY, curQuadPointParamCoord[0], curQuadPointParamCoord[1], curQuadWeight);
                            }
                            else
                            {
                                Float curQuadPointX = TransformSamplePosition(quadData1DI.nodes[k], llCorner[0], uuCorner[0]);
                                Float curQuadPointY = TransformSamplePosition(quadData1DJ.nodes[l], llCorner[1], uuCorner[1]);
                                curQuadWeight = quadData1DI.weights[k] * quadData1DJ.weights[l] * patchWeightCoeff;

                                Vec2 curQuadPoint = Vec2(curQuadPointX, curQuadPointY);
                                Vec2 suggestUV = Vec2(uMin + 0.5, vMin + 0.5);
                                curQuadPointParamCoord = this->ParametricCoordAtRest(curQuadPoint, suggestUV);

                                BSIPC_INFO("Quad Point: x: {:02.2f} y: {:02.2f} u: {:02.2f} v: {:02.2f} w: {:04.8f}",
                                    curQuadPointX, curQuadPointY, curQuadPointParamCoord[0], curQuadPointParamCoord[1], curQuadWeight);
                            }

                            // cache derivatives
                            Mat<2, 2> pd_Material_Param = this->PD_Material_Param(curQuadPointParamCoord[0], curQuadPointParamCoord[1]);
                            Mat<2, 3> pd2_Material_Param = this->PD2_Material_Param(curQuadPointParamCoord[0], curQuadPointParamCoord[1]);
                            Mat<2, 2> pd_Param_Material = this->PD_Param_Material(curQuadPointParamCoord[0], curQuadPointParamCoord[1]);
                            Mat<2, 3> pd2_Param_Material = this->PD2_Param_Material(
                                pd_Param_Material, pd_Material_Param, pd2_Material_Param
                            );

                            QuadPoint curQuad = QuadPoint
                            {
                                curQuadPointParamCoord[0], curQuadPointParamCoord[1], curQuadWeight,
                                pd_Param_Material, pd2_Param_Material, pd_Material_Param, pd2_Material_Param
                            };

                            UInt index = quadPtOffset + k * orderJ + l;

                            this->ssQuadPoints[index] = curQuad;

                            UInt curQuadIndex = index;

                            IVec2 supportStart = this->SupportStartOf(curQuad.U(), curQuad.V());
                            for (UInt n = 0; n != 9; ++n)   // for each control point whose support contains the quadrature point
                            {
                                UInt uOffset = n % 3;
                                UInt vOffset = n / 3;

                                IVec2 curSupport = supportStart + IVec2(uOffset, vOffset);
                                UInt curSupportIndex = this->GetControlPointIndex(curSupport[0], curSupport[1]);

                                BSIPC_ASSERT(
                                    curSupportIndex < this->nearbyQuadPointCache.size(),
                                    "Nearby Quad Point Cache size too small!");

                                this->nearbyQuadPointCache[curSupportIndex].insert(curQuadIndex);
                            }

                            //w += curQuadWeight;
                            //++numQuadPoints;
                        }
                }
            );
        }

        // Init bending (1-1) quadrature points
        UInt uMax = this->UDomainMax(), vMax = this->VDomainMax();

        UInt patchCnt = uMax * vMax;
        this->bdQuadPoints.resize(patchCnt, QuadPoint::Zero());

        tbb::parallel_for(
            0, static_cast<Int>(patchCnt), 1, [&](Int iter)
            //for (Int iter = 0; iter != static_cast<Int>(patchCnt); ++ iter)
            {
                UInt i = static_cast<UInt>(iter % uMax);
                UInt j = static_cast<UInt>(iter / uMax);

                Vec2 llCorner = this->EvalAtRest(i, j);
                Vec2 luCorner = this->EvalAtRest(i, j + 1);
                Vec2 ulCorner = this->EvalAtRest(i + 1, j);
                Vec2 uuCorner = this->EvalAtRest(i + 1, j + 1);

                Float patchArea = 0;
                patchArea += trigArea(llCorner, luCorner, ulCorner);
                patchArea += trigArea(luCorner, ulCorner, uuCorner);

                Float curQuadWeight = patchArea;

                Vec2 curQuadPoint = (llCorner + luCorner + ulCorner + uuCorner) / 4;
                Vec2 suggestUV = Vec2(i + 0.5, j + 0.5);
                Vec2 curQuadPointParamCoord = this->ParametricCoordAtRest(curQuadPoint, suggestUV);

                //BSIPC_INFO("Quad Point: {:02.2f} v: {:02.2f} w: {:04.4f}",
                //    curQuadPointParamCoord[0], curQuadPointParamCoord[1], curQuadWeight);

                // cache derivatives
                Mat<2, 2> pd_Material_Param = this->PD_Material_Param(curQuadPointParamCoord[0], curQuadPointParamCoord[1]);
                Mat<2, 3> pd2_Material_Param = this->PD2_Material_Param(curQuadPointParamCoord[0], curQuadPointParamCoord[1]);
                Mat<2, 2> pd_Param_Material = this->PD_Param_Material(curQuadPointParamCoord[0], curQuadPointParamCoord[1]);
                Mat<2, 3> pd2_Param_Material = this->PD2_Param_Material(
                    pd_Param_Material, pd_Material_Param, pd2_Material_Param
                );

                QuadPoint curQuad = QuadPoint
                {
                    curQuadPointParamCoord[0], curQuadPointParamCoord[1], curQuadWeight,
                    pd_Param_Material, pd2_Param_Material, pd_Material_Param, pd2_Material_Param
                };

                this->bdQuadPoints[iter] = curQuad;
            }
        );

        //BSIPC_INFO("=================================================");
        //for (auto& quad : this->quadPoints)
        //    BSIPC_INFO("Quad Point: u: {:02.2f} v: {:02.2f} w: {:02.2f}", quad.U(), quad.V(), quad.weight);

        //BSIPC_INFO("Total Weight: {:02.2f}", w);
        //BSIPC_INFO("Total Quad Points: {}", numQuadPoints);

        if (!IsQuadSchemeGlobal(scheme) && scheme != QuadScheme::LOCAL_DUAL_CHESSBOARD)
        {
            UInt totalQuadPts = this->QuadPointCnt(scheme);
            BSIPC_INFO("Computed Quad Points Cnt: {}", totalQuadPts);
        }

        //this->solver->AppendToLog("Total Quad Points: " + std::to_string(numQuadPoints) + "\n");
        //this->solver->LogSep();
    }

    UInt BSSurface::QuadPointCnt(QuadScheme scheme) const
    {
        // this function is for accelerating quadrature computation for local schemes.
        // global schemes should not appear here.

        BSIPC_ASSERT(
            !IsQuadSchemeGlobal(scheme), "Global quadrature schemes should not query about quadrature points"
        );

        if (scheme == QuadScheme::LOCAL_UNIFORM)
            return 2 * 2 * this->UDomainMax() * this->VDomainMax();
        else if (scheme == QuadScheme::LOCAL_EDGE_DENSE)
        {
            UInt edgePatchCnt = 2 * this->UDomainMax() + 2 * this->VDomainMax() - 4;
            UInt nonEdgePatchCnt = this->UDomainMax() * this->VDomainMax() - edgePatchCnt;

            // 3 * 2 * [edge/corner patch] + interior + [corner patch, they are actually 3 * 3]
            return 3 * 2 * edgePatchCnt + 2 * 2 * nonEdgePatchCnt + 3 * 4;
        }
        else if (scheme == QuadScheme::LOCAL_CHESSBOARD)
        {
            UInt edgePatchCnt = 2 * this->UDomainMax() + 2 * this->VDomainMax() - 4;
            UInt nonEdgePatchCnt = this->UDomainMax() * this->VDomainMax() - edgePatchCnt;

            return 3 * 2 * edgePatchCnt + 2 * nonEdgePatchCnt + 3 * 4;
        }
        else if (scheme == QuadScheme::LOCAL_CHESSBOARD_DENSE)
        {
            UInt edgePatchCnt = 2 * this->UDomainMax() + 2 * this->VDomainMax() - 4;
            UInt interiorQuadCnt = 0;

            UInt cols = this->UDomainMax(), rows = this->VDomainMax();
            if (rows % 2 == 0 || cols % 2 == 0)
                interiorQuadCnt = (rows - 2) * (cols - 2) / 2 * 5;
            else
                interiorQuadCnt = (cols * rows - 2 * (cols + rows) + 3) / 2 * 5 + 4;

            return 3 * 2 * edgePatchCnt + interiorQuadCnt + 3 * 4;
        }

        BSIPC_ASSERT(true, "[Fallthrough] Invalid quad scheme");
        return std::numeric_limits<UInt>::max();
    }

    UInt BSSurface::QuadPointStartIndexAtPatch(QuadScheme scheme, UInt patchIndexI, UInt patchIndexJ) const
    {
        BSIPC_ASSERT(!IsQuadSchemeGlobal(scheme), "Global quadrature schemes should not query about patch index");

        UInt uMax = this->UDomainMax();
        UInt vMax = this->VDomainMax();

        UInt globalPatchIndex = patchIndexJ * uMax + patchIndexI;

        // quad points in patches with global index <= (patchIndexI, patchIndexJ)
        UInt prevQuadCnt = tbb::parallel_reduce(
            tbb::blocked_range<UInt>(0, globalPatchIndex),
            0,
            [&](const tbb::blocked_range<UInt>& r, UInt sum) -> UInt
            {
                for (UInt idx = r.begin(); idx != r.end(); ++idx)
                {
                    UInt curIndexI = idx % uMax;
                    UInt curIndexJ = idx / uMax;

                    std::pair<UInt, UInt> localOrderRes = this->LocalPatchOrder(scheme, curIndexI, curIndexJ);

                    UInt localQuadCnt = localOrderRes.first * localOrderRes.second;
                    sum += localQuadCnt;
                }
                return sum;
            },
            [](UInt a, UInt b) { return a + b; }
        );

        return prevQuadCnt;
    }

    std::pair<UInt, UInt> BSSurface::LocalPatchOrder(QuadScheme scheme, UInt i, UInt j) const
    {
        UInt orderI = 2;
        UInt orderJ = 2;
        if (scheme == QuadScheme::LOCAL_UNIFORM)
            ;
        else if (scheme == QuadScheme::LOCAL_EDGE_DENSE)
        {
            if (i == 0 || i == this->UDomainMax() - 1)
                orderI = 3;
            if (j == 0 || j == this->VDomainMax() - 1)
                orderJ = 3;
        }
        else if (scheme == QuadScheme::LOCAL_CHESSBOARD)
        {
            if (i % 2 == j % 2)
            {
                orderI = 2;
                orderJ = 1;
            }
            else
            {
                orderI = 1;
                orderJ = 2;
            }

            // edge dense overrides chessboard
            if (i == 0 || i == this->UDomainMax() - 1)
            {
                orderI = 3;
                orderJ = (orderJ >= 2) ? orderJ : 2;
            }
            if (j == 0 || j == this->VDomainMax() - 1)
            {
                orderJ = 3;
                orderI = (orderI >= 2) ? orderI : 2;
            }
        }
        else if (scheme == QuadScheme::LOCAL_CHESSBOARD_DENSE)
        {
            // 2-by-2 surrounded by 1-1
            if (i % 2 == j % 2)
            {
                orderI = 2;
                orderJ = 2;
            }
            else
            {
                orderI = 1;
                orderJ = 1;
            }

            // edge dense overrides chessboard
            if (i == 0 || i == this->UDomainMax() - 1)
            {
                orderI = 3;
                orderJ = (orderJ >= 2) ? orderJ : 2;
            }
            if (j == 0 || j == this->VDomainMax() - 1)
            {
                orderJ = 3;
                orderI = (orderI >= 2) ? orderI : 2;
            }
        }

        return std::make_pair(orderI, orderJ);
    }

    QuadNodesData BSSurface::TrimDenseQuadPoints(const std::map<UInt, UInt>& specLimits, const QuadNodesData& data, UInt limit, UInt cellCnt) const
    {
        QuadNodesData trimmedData;
        UInt tunedLimit = limit;
        UInt quadOrder = static_cast<UInt>(data.nodes.size());

        for (UInt intervalStart = 0; intervalStart != cellCnt; ++intervalStart)
        {
            UInt lowerBoundIndex = 0, upperBoundIndex = 0;
            Bool lowerBoundFound = false, upperBoundFound = false;

            // Identify nodes with value in [intervalStart, intervalStart + 1)
            for (UInt i = 0; i != data.nodes.size(); ++i)
            {
                if (data.nodes[i] >= intervalStart && !lowerBoundFound)
                {
                    lowerBoundIndex = i;
                    lowerBoundFound = true;
                }

                if (data.nodes[i] > intervalStart + 1 && !upperBoundFound)
                {
                    upperBoundIndex = i - 1;
                    upperBoundFound = true;
                    break;
                }
            }

            // If no node exists in the interval
            if (upperBoundIndex - lowerBoundIndex < 0)
                continue;

            // If intervalStart == domainMax - 1
            if (!upperBoundFound)
                upperBoundIndex = quadOrder - 1;

            BSIPC_ASSERT(upperBoundIndex < data.nodes.size(), "Upper bound index out of bound");

            std::vector<Float> curSegmentNodes, curSegmentWeights;
            curSegmentNodes = std::vector<Float>(data.nodes.begin() + lowerBoundIndex, data.nodes.begin() + upperBoundIndex + 1);
            curSegmentWeights = std::vector<Float>(data.weights.begin() + lowerBoundIndex, data.weights.begin() + upperBoundIndex + 1);

            UInt bound = upperBoundIndex - lowerBoundIndex;

            if (specLimits.contains(intervalStart))
                tunedLimit = specLimits.at(intervalStart);
            else
                tunedLimit = limit;

            // Only merge for segments having at least two node
            if (curSegmentNodes.size() >= 2)
            {
                for (UInt iter = tunedLimit; iter <= bound; ++iter)
                {
                    // Identify the closest nodes
                    Float smallestDiff = std::numeric_limits<Float>::max();
                    UInt smallestDiffIndex = 0;

                    for (UInt curIndex = 0; curIndex != curSegmentNodes.size() - 1; ++curIndex)
                    {
                        Float diff = curSegmentNodes[curIndex + 1] - curSegmentNodes[curIndex];
                        if (diff < smallestDiff)
                        {
                            smallestDiff = diff;
                            smallestDiffIndex = curIndex;
                        }
                    }

                    // Merge them
                    MergeConsecutiveNodes(curSegmentNodes, curSegmentWeights, smallestDiffIndex);
                }
            }

            for (Float node : curSegmentNodes)
                trimmedData.nodes.push_back(node);
            for (Float weight : curSegmentWeights)
                trimmedData.weights.push_back(weight);
        }
        BSIPC_PEEK(trimmedData.nodes.size());

        return trimmedData;
    }

    void BSSurface::UpdateQuadPointCache()
    {
        Int size = static_cast<Int>(this->ssQuadPoints.size());

        tbb::parallel_for(
            0, size, 1,
            [&](UInt i)
            {
                QuadPoint& quadPoint = this->ssQuadPoints[i];
                quadPoint.deformationGradient = this->DeformationGradientAt(quadPoint.U(), quadPoint.V(), quadPoint.PD_Param_Material());
                quadPoint.pd_World_Param = this->PD_World_Param(quadPoint.U(), quadPoint.V());
                quadPoint.pd2_World_Param = this->PD2_World_Param(quadPoint.U(), quadPoint.V());
            }
        );
    }

    void BSSurface::InitMassMatrix(Float density)
    {
        const std::vector<QuadPoint>& ssQuadPoints = this->GetSsQuadPoints();
        UInt size = this->GetControlPointCnt();
        BSIPC_PEEK(size);
        this->massMatDiagEntries.resize(size, 0);

        const auto& nearbyCache = this->GetNearbyQuadPointCache();

        // TODO possibly could compute only half entries... to be optimized if needed.
        tbb::parallel_for(
            0, static_cast<Int>(size * size), 1, [&](Int i)
            //for (Int i = 0; i != size * size; ++i)
            {
                UInt matRowIndex = i % size;
                UInt matColIndex = i / size;

                Float curMassMatEntry = 0.;
                const std::set<UInt>& rowNearbyIndices = nearbyCache[matRowIndex];
                const std::set<UInt>& colNearbyIndices = nearbyCache[matColIndex];

                std::set<UInt> nearbyIndices = rowNearbyIndices;
                nearbyIndices.insert(colNearbyIndices.begin(), colNearbyIndices.end());

                for (UInt quadIndex : nearbyIndices)
                {
                    const QuadPoint& quadPoint = ssQuadPoints[quadIndex];
                    Float addComp = density * quadPoint.Weight() *
                        this->EvalCoeffAt(matRowIndex, quadPoint.U(), quadPoint.V()) *
                        this->EvalCoeffAt(matColIndex, quadPoint.U(), quadPoint.V());
                    //BSIPC_INFO("ri: {}, ci: {} qindex: {}, addComp: {}", matRowIndex, matColIndex, quadIndex, addComp);
                    curMassMatEntry += addComp;
                }

                this->massMatDiagEntries[matRowIndex] += curMassMatEntry;
            }
        );

        // {TEST} the following two terms should match
        Float massMatEntriesSum = 0.;
        for (UInt i = 0; i != size; ++i)
            massMatEntriesSum += this->massMatDiagEntries[i];

        Float massOverSurface = density * this->surfaceArea;

        BSIPC_INFO("Sum of Mass Matrix Entries: {}", massMatEntriesSum);
        BSIPC_INFO("Surface mass (density * area): {}", massOverSurface);
        
        Float error = std::abs(massMatEntriesSum - massOverSurface) / massOverSurface;
        BSIPC_INFO("mass relative error: {:.6e}", error);
        
    }

    void BSSurface::WriteQuadCache(String path) const
    {
        fkyaml::node root(fkyaml::node_type::MAPPING);

        std::vector<fkyaml::node> ssQuadPointsNode;
        for (const QuadPoint& quadPoint : this->ssQuadPoints)
            ssQuadPointsNode.push_back(fkyaml::node{ quadPoint.U(), quadPoint.V(), quadPoint.Weight() });

        std::vector<fkyaml::node> bdQuadPointsNode;
        for (const QuadPoint& quadPoint : this->bdQuadPoints)
            bdQuadPointsNode.push_back(fkyaml::node{ quadPoint.U(), quadPoint.V(), quadPoint.Weight() });

        std::vector<fkyaml::node> patchStartIndices;
        for (UInt i = 0; i != this->GetPatchStartIndices().size(); ++i)
            patchStartIndices.push_back(this->GetPatchStartIndices()[i]);

        std::vector<fkyaml::node> massMatrixEntries;
        for (UInt i = 0; i != this->GetMassMatEntries().size(); ++i)
            massMatrixEntries.push_back(this->GetMassMatEntries()[i]);

        root[Config::SS_QUAD_POINTS_STR] = ssQuadPointsNode;
        root[Config::BD_QUAD_POINTS_STR] = bdQuadPointsNode;
        root[Config::PATCH_START_INDEX_STR] = patchStartIndices;
        root[Config::MASS_MATRIX_STR] = massMatrixEntries;

        std::ofstream ofs(path);
        if (ofs)
            ofs << root;
    }

    void BSSurface::ReadQuadData(String path)
    {
        std::ifstream ifs(path);

        try
        {
            if (!ifs)
            {
                std::string msg = "error opening config file " + path;
                throw fkyaml::exception(msg.c_str());
            }
            fkyaml::node root = fkyaml::node::deserialize(ifs);

            auto ssQuadPointsNode = root;
            LoadNodeFromYaml(ssQuadPointsNode, root, Config::SS_QUAD_POINTS_STR);

            this->ssQuadPoints.clear();
            for (const fkyaml::node& node : ssQuadPointsNode)
            {
                BSIPC_ASSERT(node.is_sequence() && node.size() == 3, "Quad point node is not a sequence or has size different from 3");

                Float u = node[0].get_value<Float>(), v = node[1].get_value<Float>(), weight = node[2].get_value<Float>();

                Mat<2, 2> pd_Material_Param = this->PD_Material_Param(u, v);
                Mat<2, 3> pd2_Material_Param = this->PD2_Material_Param(u, v);
                Mat<2, 2> pd_Param_Material = this->PD_Param_Material(u, v);
                Mat<2, 3> pd2_Param_Material = this->PD2_Param_Material(
                    pd_Param_Material, pd_Material_Param, pd2_Material_Param
                );

                QuadPoint curQuadPoint = QuadPoint(u, v, weight, pd_Param_Material, pd2_Param_Material, pd_Material_Param, pd2_Material_Param);
                this->ssQuadPoints.push_back(curQuadPoint);
            }

            auto bdQuadPointsNode = root;
            LoadNodeFromYaml(bdQuadPointsNode, root, Config::BD_QUAD_POINTS_STR);

            this->bdQuadPoints.clear();
            for (const fkyaml::node& node : bdQuadPointsNode)
            {
                BSIPC_ASSERT(node.is_sequence() && node.size() == 3, "Quad point node is not a sequence or has size different from 3");

                Float u = node[0].get_value<Float>(), v = node[1].get_value<Float>(), weight = node[2].get_value<Float>();

                Mat<2, 2> pd_Material_Param = this->PD_Material_Param(u, v);
                Mat<2, 3> pd2_Material_Param = this->PD2_Material_Param(u, v);
                Mat<2, 2> pd_Param_Material = this->PD_Param_Material(u, v);
                Mat<2, 3> pd2_Param_Material = this->PD2_Param_Material(
                    pd_Param_Material, pd_Material_Param, pd2_Material_Param
                );

                QuadPoint curQuadPoint = QuadPoint(u, v, weight, pd_Param_Material, pd2_Param_Material, pd_Material_Param, pd2_Material_Param);
                this->bdQuadPoints.push_back(curQuadPoint);
            }

            auto patchStartIndicesNode = root;
            LoadNodeFromYaml(patchStartIndicesNode, root, Config::PATCH_START_INDEX_STR);

            for (const fkyaml::node& node : patchStartIndicesNode)
            {
                UInt patchStartIndex = node.get_value<UInt>();
                this->patchQuadStartIndices.push_back(patchStartIndex);
            }

            auto massMatrixNode = root;
            LoadNodeFromYaml(massMatrixNode, root, Config::MASS_MATRIX_STR);

            for (const fkyaml::node& node : massMatrixNode)
            {
                Float massMatrixEntry = node.get_value<Float>();
                this->massMatDiagEntries.push_back(massMatrixEntry);
            }
        }
        catch (fkyaml::exception& e)
        {
            BSIPC_WARN(e.what());
            BSIPC_WARN("Error loading quadrature point cache");
        }
    }

    void BSSurface::TestQuadData(String path)
    {
        std::ifstream ifs(path);

        try
        {
            if (!ifs)
            {
                std::string msg = "error opening config file " + path;
                throw fkyaml::exception(msg.c_str());
            }
            fkyaml::node root = fkyaml::node::deserialize(ifs);

            auto ssQuadPointsNode = root;
            LoadNodeFromYaml(ssQuadPointsNode, root, Config::SS_QUAD_POINTS_STR);

            std::vector<QuadPoint> ssTestQuadPoints;
            ssTestQuadPoints.clear();
            for (const fkyaml::node& node : ssQuadPointsNode)
            {
                BSIPC_ASSERT(node.is_sequence() && node.size() == 3, "Quad point node is not a sequence or has size different from 3");

                Float u = node[0].get_value<Float>(), v = node[1].get_value<Float>(), weight = node[2].get_value<Float>();

                Mat<2, 2> pd_Material_Param = this->PD_Material_Param(u, v);
                Mat<2, 3> pd2_Material_Param = this->PD2_Material_Param(u, v);
                Mat<2, 2> pd_Param_Material = this->PD_Param_Material(u, v);
                Mat<2, 3> pd2_Param_Material = this->PD2_Param_Material(
                    pd_Param_Material, pd_Material_Param, pd2_Material_Param
                );

                QuadPoint curQuadPoint = QuadPoint(u, v, weight, pd_Param_Material, pd2_Param_Material, pd_Material_Param, pd2_Material_Param);
                ssTestQuadPoints.push_back(curQuadPoint);
            }

            auto bdQuadPointsNode = root;
            LoadNodeFromYaml(bdQuadPointsNode, root, Config::BD_QUAD_POINTS_STR);

            std::vector<QuadPoint> bdTestQuadPoints;
            bdTestQuadPoints.clear();
            for (const fkyaml::node& node : bdQuadPointsNode)
            {
                BSIPC_ASSERT(node.is_sequence() && node.size() == 3, "Quad point node is not a sequence or has size different from 3");

                Float u = node[0].get_value<Float>(), v = node[1].get_value<Float>(), weight = node[2].get_value<Float>();

                Mat<2, 2> pd_Material_Param = this->PD_Material_Param(u, v);
                Mat<2, 3> pd2_Material_Param = this->PD2_Material_Param(u, v);
                Mat<2, 2> pd_Param_Material = this->PD_Param_Material(u, v);
                Mat<2, 3> pd2_Param_Material = this->PD2_Param_Material(
                    pd_Param_Material, pd_Material_Param, pd2_Material_Param
                );

                QuadPoint curQuadPoint = QuadPoint(u, v, weight, pd_Param_Material, pd2_Param_Material, pd_Material_Param, pd2_Material_Param);
                bdTestQuadPoints.push_back(curQuadPoint);
            }

            std::vector<UInt> patchQuadStartIndicesTest;
            std::vector<Float> massMatDiagEntriesTest;

            auto patchStartIndicesNode = root;
            LoadNodeFromYaml(patchStartIndicesNode, root, Config::PATCH_START_INDEX_STR);

            for (const fkyaml::node& node : patchStartIndicesNode)
            {
                UInt patchStartIndex = node.get_value<UInt>();
                patchQuadStartIndicesTest.push_back(patchStartIndex);
            }

            auto massMatrixNode = root;
            LoadNodeFromYaml(massMatrixNode, root, Config::MASS_MATRIX_STR);

            for (const fkyaml::node& node : massMatrixNode)
            {
                Float massMatrixEntry = node.get_value<Float>();
                massMatDiagEntriesTest.push_back(massMatrixEntry);
            }

            if (ssTestQuadPoints != this->ssQuadPoints)
            {
                if (ssTestQuadPoints.size() != this->ssQuadPoints.size())
                    BSIPC_WARN("size of ssQuad unequal");

                for (UInt i = 0; i != ssTestQuadPoints.size(); ++i)
                {
                    if (ssTestQuadPoints[i] != this->ssQuadPoints[i])
                    {
                        BSIPC_INFO("ssQuadPoint[{}] unequal", i);
                        BSIPC_WARN("ssTestQuadPoint: {}", ssTestQuadPoints[i].Info());
                        BSIPC_WARN("ssQuadPoint: {}", this->ssQuadPoints[i].Info());
                    }
                }
            }
            
            if (bdTestQuadPoints != this->bdQuadPoints)
                BSIPC_WARN("bdQuadPoints are not equal");

            if (patchQuadStartIndicesTest != this->patchQuadStartIndices)
                BSIPC_WARN("patchQuadStartIndices are not equal");

            if (massMatDiagEntriesTest != this->massMatDiagEntries)
                BSIPC_WARN("massMatDiagEntries are not equal");
        }
        catch (fkyaml::exception& e)
        {
            BSIPC_WARN(e.what());
            BSIPC_WARN("Error loading quadrature point cache");
        }
    }
}
