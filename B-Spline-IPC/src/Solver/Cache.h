#pragma once
#include "bsfwd.h"

namespace BSIPC
{
    // Stores variable that can be reused in each step, between different functions
    // all fields should be defined as std::optional, and Reset() should be called at the beginning of each step
    struct StepCache
    {
    public:
        StepCache() = default;
        StepCache(UInt startStep) : stepCnt(startStep) { }

        inline UInt CurStep() const { return this->stepCnt; }
        inline void StartFromStep(UInt step) 
        { 
            if (this->stepCnt != 0)
                BSIPC_CRITICAL("Attempting to modify step number after simulation has launched!");
            this->stepCnt = step; 
        }

        inline void NextStep(UInt velSize)
        {
            estVels.clear();
            estVels.resize(velSize, Vec3::Zero());

            ++this->stepCnt;
            this->newtonIterCnt = 0;
        }

    public:
        UInt newtonIterCnt = 0;
        std::vector<Vec3> estVels;
        UInt lineSearchIterCnt = 0;

    private:
        UInt stepCnt = 0;
    };

    // Stores variable that can be reused in each quadrature point, but can vary in each step
    // all fields should be defined as std::optional, and Reset() should be called at each change of quadrature point
    class QuadratureCache
    {
    public:
        QuadratureCache() = default;

        // NOTE cannot update F and J here as they depend on the target B-Spline surface, which is only linked to the solver
        Mat<3, 2> F;                        // Deformation Gradient
        Mat<2, 2> PD_Param_Coord;           // Partial Derivative of Parametric Coordinates
        Float J = 0.f;                      // Determinant of Deformation Gradient

        /**
         * @params (u, v) are the next quadrature point for which the cache is to be used
         */
        inline void NextQuadPoint(Float u, Float v)
        {
            this->F = Mat<3, 2>::Zero();
            this->J = 0.;
            pd_dgCtrlpt.clear();

            quadPtIndex = std::make_pair(u, v);
        }

    private:
        std::pair<Float, Float> quadPtIndex;
        std::unordered_map<uint64_t, Vec2> pd_dgCtrlpt;             // Vec2 = (Gamma^1_ij, Gamma^2_ij)
    };

    struct AnimatedMeshCache
    {
        UInt startIdx = std::numeric_limits<UInt>::max();
        UInt vertCnt = std::numeric_limits<UInt>::max();
        std::vector<Vec3> curPos;
    };

}
