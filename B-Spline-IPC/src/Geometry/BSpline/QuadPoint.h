#pragma once
#include "bsfwd.h"

#include "Utility/EigenHelper.h"

namespace BSIPC
{
    struct QuadPoint
    {
    private:
        Float u = 0;
        Float v = 0;
        Float weight = 0;
        
        // Quantities that are related only to rest shape are marked as const, and should be computed during initialization

        // Material space related PDs, should be initialized when initializing the quad points
        Mat<2, 2> pd_Param_Material             = Mat<2, 2>::Zero();
        Mat<2, 3> pd2_Param_Material            = Mat<2, 3>::Zero();

        Mat<2, 2> pd_Material_Param             = Mat<2, 2>::Zero();
        Mat<2, 3> pd2_Material_Param            = Mat<2, 3>::Zero();

    public:
        // World space related PDs, updated each time the control points are updated
        Mat<3, 2> pd_World_Param                = Mat<3, 2>::Zero();
        Mat<3, 3> pd2_World_Param               = Mat<3, 3>::Zero();

        Mat<3, 2> deformationGradient           = Mat<3, 2>::Zero();

        QuadPoint(
            Float u, Float v, Float weight,
            const Mat<2, 2>& pd_Param_Material,
            const Mat<2, 3>& pd2_Param_Material,
            const Mat<2, 2>& pd_Material_Param,
            const Mat<2, 3>& pd2_Material_Param
        ) : 
            u(u), v(v), weight(weight),
            pd_Param_Material(pd_Param_Material),
            pd2_Param_Material(pd2_Param_Material),
            pd_Material_Param(pd_Material_Param),
            pd2_Material_Param(pd2_Material_Param)
        {}

        QuadPoint& operator=(const QuadPoint& other)
        {
            this->u = other.U();
            this->v = other.V();
            this->weight = other.Weight();
            this->pd_Param_Material = other.PD_Param_Material();
            this->pd2_Param_Material = other.PD2_Param_Material();
            this->pd_Material_Param = other.PD_Material_Param();
            this->pd2_Material_Param = other.PD2_Material_Param();
            this->pd_World_Param = other.pd_World_Param;
            this->pd2_World_Param = other.pd2_World_Param;
            this->deformationGradient = other.deformationGradient;

            return *this;
        }

        inline String Info() const
        {
            return fmt::format("== QuadPoint ({}, {}) with weight {} ==\n\
                pd_P_M: {}\npd2_P_M: {}\npd_M_P: {}\npd2_M_P: {}\npd_W_P: {}\npd2_W_P: {}\nF: {}",
                u, v, weight, ToStr(pd_Param_Material), ToStr(pd2_Param_Material), ToStr(pd_Material_Param), ToStr(pd2_Material_Param), ToStr(pd_World_Param), ToStr(pd2_World_Param), ToStr(deformationGradient));
        }

        Float U() const { return u; }
        Float V() const { return v; }
        Float Weight() const { return weight; }
        const Mat<2, 2>& PD_Param_Material() const { return pd_Param_Material; }
        const Mat<2, 3>& PD2_Param_Material() const { return pd2_Param_Material; }
        const Mat<2, 2>& PD_Material_Param() const { return pd_Material_Param; }
        const Mat<2, 3>& PD2_Material_Param() const { return pd2_Material_Param; }

        static inline QuadPoint Zero()
        {
            return QuadPoint(0, 0, 0, Mat<2, 2>::Zero(), Mat<2, 3>::Zero(), Mat<2, 2>::Zero(), Mat<2, 3>::Zero());
        }
    };

    inline Bool operator==(const QuadPoint& lhs, const QuadPoint& rhs)
    {
        return (lhs.U() == rhs.U() && lhs.V() == rhs.V() && lhs.Weight() == rhs.Weight()
            && lhs.PD_Param_Material() == rhs.PD_Param_Material()
            && lhs.PD2_Param_Material() == rhs.PD2_Param_Material()
            && lhs.PD_Material_Param() == rhs.PD_Material_Param()
            && lhs.PD2_Material_Param() == rhs.PD2_Material_Param()
            && lhs.pd_World_Param == rhs.pd_World_Param
            && lhs.pd2_World_Param == rhs.pd2_World_Param
            && lhs.deformationGradient == rhs.deformationGradient);
    }
}
