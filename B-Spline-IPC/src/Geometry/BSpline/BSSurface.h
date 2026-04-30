#pragma once

#include "bsfwd.h"

#include "QuadPoint.h"
#include "Utility/Math/GaussianQuadrature.h"
#include "Solver/SolverConfig.h"

namespace BSIPC
{
    class Solver;

    class BSSurface
    {
        /* ***************************** BS Surface Layout *****************************
         * The nodes are indexed with
         * (0, 0) -- ... -- (0, w - 1) --
         * (1, 0) -- ... -- (1. w - 1) -- 
         * ... / ... -- 
         * (h - 1, 0) -- ... -- (h - 1, w - 1)
         * */

        // This results in [[0, w-1]] x [[0, h-1]]

    public:
        BSSurface(SolverConfig config, UInt index, const BSSurfaceInfo& info);

        /* ***************************** BS Surface Evaluations ***************************** */

        /**** Parameteric Evaluation of B-Spline Surfaces ****
         * B_2^2(u, v) = \sum_ij N_i(u) N_j(v) P_{i, j}
         */

        // Computes \sum_ij N_i(u) N_j(v) P_{i, j}
        Vec3 EvalAt(Float u, Float v) const;

        // Computes \sum_ij N_i(u) N_j(v) P_{i, j}^0 where P_{**}^0 is the rest configuration
        Vec2 EvalAtRest(Float u, Float v) const;

        // Computes N_i(u) N_j(v) which is the coefficient of evaluation at (u, v) w.r.t. (i, j)-th control point
        Float EvalCoeffAt(UInt i, UInt j, Float u, Float v) const;
        // Overload, (i, j) corresponds to (is replaced by) the single dimensional index
        Float EvalCoeffAt(UInt index, Float u, Float v) const;

        // Computes the normal at (u, v).
        Vec3 NormalAt(Float u, Float v) const;

        // Computes the distance between the point whose parametric coordinates is (u, v) and the point at pos.
        Float EvalDistancebetween(Float u, Float v, Vec2 pos) const;

        // TODO modify this if the element is not rectangle
        inline Float LengthSideX() const;
        inline Float LengthSideY() const;

        /**
         * Computes the partial derivative dr/du and dr/dv where r is the distance between (u, v) in parametric coordinates and pos in material space.
         * ( dr/du dr/dv )
         * */
        Vec2 PD_Distance_ParamCoord(Float u, Float v, Vec2 pos) const;

        // Computes the parametric coordinate (u, v) at the point X in the material space, at rest.
        Vec2 ParametricCoordAtRest(Vec2 pos, Vec2 suggestUV) const;

        /* ***************************** Partial Derivatives ***************************** */

        Mat<3, 2> DeformationGradientAt(Float u, Float v, const Mat<2, 2>& pd_Param_Material);

        //// Computes the function for Gamma^{n}_{i, j} = partial derivative of F_mn w.r.t. P_ij at (u, v).
        //// m does not need to passed in here as the derivative is the same for every point, independent of which entry in the column
        //// For convenience of caching in solver, returns (Gamma^1_ij, Gamma^2_ij)
        //Vec2 PDDefgradCtrlpt(UInt i, UInt j, Float u, Float v) const;

        /**
         * Computes the partial derivative of the parametric coordinates w.r.t. the material space coordinates.
         *
         * For coordinates (X1, X2) in material space and (u, v) in parametric space, returns the matrix:
         * ( du/dX1  du/dX2 )
         * ( dv/dX1  dv/dX2 )
         */
        Mat<2, 2> PD_Param_Material(Float u, Float v) const;

        /**
         * Computes the second order partial derivative of parametric coordinates w.r.t, the material space coordinates.
         * REQUIRES all other derivatives in the quadpoints cache have been updated
         *          computing the second order partial derivatives will use other derivatives
         * 
         * For coordinates (X1, X2) in material space and (u, v) in parametric space, returns the matrix:
         * ( d^2u/dX1^2  d^2u/dX1dX2  d^2u/dX2^2 )
         * ( d^2v/dX1^2  d^2v/dX1dX2  d^2v/dX2^2 )
         */
        Mat<2, 3> PD2_Param_Material(Mat<2, 2> PD_Param_Material, Mat<2, 2> PD_Material_Param, Mat<2, 3> PD2_Material_Param) const;

        /**
         * Computes the partial derivative of the material space coordinates w.r.t. the parametric coordinates.
         *
         * For coordinates (X1, X2) in material space and (u, v) in parametric space, returns the matrix:
         * ( dX1/du  dX2/du )
         * ( dX1/dv  dX2/dv )
         */
        Mat<2, 2> PD_Material_Param(Float u, Float v) const;

        /**
         * Computes the second order partial derivative of the material space coordinates w.r.t. the parametric coordinates.
         * 
         * For coordinates (X1, X2) in material space and (u, v) in parametric space, returns the matrix:
         * ( d^2X1/du^2  d^2X1/duv  d^2X1/dv^2 )
         * ( d^2X2/du^2  d^2X2/duv  d^2X2/dv^2 )
         */
        Mat<2, 3> PD2_Material_Param(Float u, Float v) const;

        /**
         * Computes the partial derivative of the world space coordinates w.r.t. the parametric coordinates.
         * 
         * For coordinates (x, y, z) in the world space and (u, v) in parametric space, returns the matrix:
         * ( dx/du  dx/dv )
         * ( dy/du  dy/dv )
         * ( dz/du  dz/dv )
         */
        Mat<3, 2> PD_World_Param(Float u, Float v) const;

        /**
         * Computes the second order partial derivative of the world space coordinates w.r.t. the parametric coordinates.
         * 
         * For coordinates (x, y, z) in the world space and (u, v) in parametric space, returns the matrix:
         * ( d^2x/du^2  d^2x/duv  d^2x/dv^2 )
         * ( d^2y/du^2  d^2y/duv  d^2y/dv^2 )
         * ( d^2z/du^2  d^2z/duv  d^2z/dv^2 )
         */
        Mat<3, 3> PD2_World_Param(Float u, Float v) const;

        /**
         * dF/dP_ij for a fixed P_ij. The result should be a 3-by-6 matrix, but since the rows are identical after being shifted 
         *    it is sufficient to return the row vector.
         * The result should look like
         * 
         * | dF11/dPx dF12/dPx dF21/dPx dF22/dPx dF31/dPx dF32/dPx |
         * | dF11/dPy dF12/dPy dF21/dPy dF22/dPy dF31/dPy dF32/dPy |
         * | dF11/dPz dF12/dPz dF21/dPz dF22/dPz dF31/dPz dF32/dPz |
         * =
         * | dF11/dPx dF12/dPx                                     |
         * |                   dF21/dPy dF22/dPy                   | 
         * |                                     dF31/dPz dF32/dPz |
         * where dF11/dPx = dF21/dPy = dF31/dPz, and same for column 2
         * 
         * Returns Vec2(dF11/dPx, dF12/dPx)
         */
        Vec2 PD_DeformGrad_CtrlPt_At(UInt i, UInt j, Float u, Float v) const;

        /* ***************************** Quadrature Preparation ***************************** */

        // Precompute parametric coordinates of quadrature points based on solver configuration
        // Precompute also all material space related PDs
        void PrepareQuadPoints();

        // Query indices for parallelizing computation of local quadrature points
        UInt QuadPointCnt(QuadScheme scheme) const;
        UInt QuadPointStartIndexAtPatch(QuadScheme scheme, UInt patchIndexI, UInt patchIndexJ) const;
        std::pair<UInt, UInt> LocalPatchOrder(QuadScheme scheme, UInt patchIndexI, UInt patchIndexJ) const;

        /** 
         * Computes the trimmed (merged) quadrature points according to the specification
         * 
         * @params specLimits is a map from the cell index to maximum number of quad points in that cell
         * @params data is the precomputed quad points, without any trimming. It should be in world coordinates (material space)
         *         after transformation
         * @params limit is the maximum number of quad points in the cell, without being specified in specLimits
         * @params cellCnt is the total number of cells in the surface, along a certain direction
         */
        QuadNodesData TrimDenseQuadPoints(const std::map<UInt, UInt>& specLimits, const QuadNodesData& data, UInt limit, UInt cellCnt) const;

        // Update world space related PDs in the quad points
        void UpdateQuadPointCache();

        /* ***************************** Init Mass Matrix ***************************** */

        void InitMassMatrix(Float density);

        /* ***************************** Coordinate Indices ***************************** */
         
        inline IVec2 GetControlPointPosIndex(UInt index) const
        { 
            BSIPC_ASSERT(index < this->controlPoints.size(), "Index out of bounds")
            return IVec2(index % this->resolutionX, index / this->resolutionX);
        }

        inline UInt GetControlPointIndex(UInt i, UInt j) const 
        { 
            //BSIPC_INFO("i: {}, j: {}, w: {}, h: {}", i, j, this->resolutionX, this->resolutionY);
            BSIPC_ASSERT(i < this->resolutionX && j < this->resolutionY, "Index out of bounds")
            return j * this->resolutionX + i;
        }

        // returns the index of the starting point (top-left, minimum of u-/v-axis)
        inline IVec2 GetPatchPosIndex(UInt index) const
        {
            BSIPC_ASSERT(index < this->GetPatchCnt(), "Index out of bounds")
            return IVec2(index % this->UDomainMax(), index / this->UDomainMax());
        }

        inline UInt GetPatchIndex(UInt i, UInt j) const
        {
            BSIPC_ASSERT(i < this->UDomainMax() && j < this->VDomainMax(), "Index out of bounds")
            return j * this->UDomainMax() + i;
        }

#if defined BSIPC_NUMERIC_TEST
        inline Vec3& GetControlPoint(UInt index) { return this->controlPoints[index]; }
#endif

        inline Vec3 GetControlPoint(UInt i, UInt j) const { return this->controlPoints[j * this->resolutionX + i]; }     // i-th row, j-th column
        inline Vec2 GetRestControlPoint(UInt i, UInt j) const { return this->restControlPoints[j * this->resolutionX + i]; }

        /* ***************************** Support ***************************** */

        // Returns the support of the basis function N_i(u) for the i-th control point.
        // `i` is the single-dimensional index, not the global index of control point on the surface
        inline Vec2 SingleDimSupport(UInt i, UInt bound) const;

        /**
         * @params i, j are the indices of the control point
         * @params u, v are the parametric coordinates
         *
         * Returns true if the corresponding basis function evaluates to nonzero value on (u, v).
         */
        inline Bool CoordInNodeOpenSupport(UInt i, UInt j, Float u, Float v) const;
        inline Bool CoordInNodeClosedSupport(UInt i, UInt j, Float u, Float v) const;

        // Given (u, v), returns the set of control points whose corresponding bases have nonzero evaluation at (u, v)
        // Returns an IVec2, corresponding to the starting index of u and v
        // Notice that by local support, there must be 3 nodes whose support contains (u, v), in both directions
        inline IVec2 SupportStartOf(Float u, Float v) const 
        { 
            return IVec2(
                std::min(static_cast<Int>(std::floor(u)), static_cast<Int>(this->UDomainMax() - 1)), 
                std::min(static_cast<Int>(std::floor(v)), static_cast<Int>(this->VDomainMax() - 1))
            ); 
        }
        inline std::vector<UInt> SupportedNodesAt(Float u, Float v) const;
        inline std::vector<UInt> SupportedNodesAt(const QuadPoint& quadPoint) const { return this->SupportedNodesAt(quadPoint.U(), quadPoint.V()); }

        std::array<BSBasisFn, 3> SelectSupportedBasisFn(Float u, UInt bound, std::array<BSBasisFn, 5> basis) const;  // Using parametric coordinates
        std::array<UInt, 3> SelectSupportedOffset(Float u) const;                                                    // Using parametric coordinates

        // Using control point index (single dimensional)
        // `bound` gives the maximum index along the parametric coordinate of consideration
        inline BSBasisFn SelectBasisFnByIndex(UInt i, const std::array<BSBasisFn, 5>& basis, UInt bound) const;                

        /* ***************************** Misc Getter/Setters ***************************** */

        inline UInt GetResolutionX() const { return this->resolutionX; }
        inline UInt GetResolutionY() const { return this->resolutionY; }
        inline UInt GetLinearResolutionX() const { return this->linearResolutionX; }
        inline UInt GetLinearResolutionY() const { return this->linearResolutionY; }
        inline UInt UDomainMax() const { return this->resolutionX - 2; }
        inline UInt VDomainMax() const { return this->resolutionY - 2; }

        inline const std::vector<QuadPoint>& GetSsQuadPoints() const { return this->ssQuadPoints; }
        inline const std::vector<QuadPoint>& GetBdQuadPoints() const { return this->bdQuadPoints; }
        inline const std::vector<UInt>& GetPatchStartIndices() const { return this->patchQuadStartIndices; }

        inline const std::vector<Vec3>& GetControlPoints() const { return this->controlPoints; }
        inline const std::vector<Vec2>& GetRestControlPoints() const { return this->restControlPoints; }
        inline UInt GetControlPointCnt() const { return static_cast<UInt>(this->controlPoints.size()); }
        inline UInt GetPatchCnt() const { return static_cast<UInt>(this->UDomainMax() * this->VDomainMax()); }

        inline Float GetSurfaceArea() const { return this->surfaceArea; }

        inline const std::vector<Float>& GetMassMatEntries() const { return this->massMatDiagEntries; }
        
        inline void UpdateControlPoints(const std::vector<Vec3>& updatedControlPoints)
        {
            BSIPC_ASSERT(this->controlPoints.size() == updatedControlPoints.size(), "Update of control points should not alter its size");
            this->controlPoints = updatedControlPoints; 
            this->UpdateQuadPointCache();
        }

        inline void BindSolver(const Solver* solver) { this->solver = solver; }

        /* ***************************** Yaml Logging Utils ***************************** */

        void WriteQuadCache(String path) const;
        void ReadQuadData(String path);
        void TestQuadData(String path);

    private:
        inline const std::vector<std::set<UInt>>& GetNearbyQuadPointCache() const { return this->nearbyQuadPointCache; }

    private:
        const Solver* solver;

        // i-th entry of the outer vector records the indices of the control points that are in the support of the i-th quadrature point
        //  initialized (resized) in InitSimulation()
        // Updated in PrepareQuadPoints()
        std::vector<std::set<UInt>> nearbyQuadPointCache;
        
        Bool irregular;
        UInt resolutionX, resolutionY;
        UInt linearResolutionX, linearResolutionY;
        Float surfaceArea;
        std::vector<Vec2> restControlPoints;                    // Control Points in the Material Space
        std::vector<Vec3> controlPoints;                        // Control Points in the World Space, subject to outside modifications
        std::vector<QuadPoint> ssQuadPoints;                    // Quadrature points for stretching, with higher resolution
        std::vector<QuadPoint> bdQuadPoints;                    // Quadrature points for bending, 1-1 per patch
        std::vector<UInt> patchQuadStartIndices;                // map: [row major index of patch] -> [start index of quad points in this patch]
        std::vector<Float> massMatDiagEntries;                  // Diagonal entries of the mass matrix
    }; 

    inline std::vector<UInt> BSSurface::SupportedNodesAt(Float u, Float v) const
    {
        std::vector<UInt> indices{};
        IVec2 supportStartIndices = this->SupportStartOf(u, v);

        //BSIPC_WARN("u: {}, v: {}, start: ({}, {})", u, v, supportStartIndices[0], supportStartIndices[1]);

        BSIPC_ASSERT(
            supportStartIndices[0] >= 0 && supportStartIndices[1] >= 0,
            "Support starting indices are negative"
        );

        for (UInt vInd = static_cast<UInt>(supportStartIndices[1]); vInd != static_cast<UInt>(supportStartIndices[1]) + 3; ++vInd)
            for (UInt uInd = static_cast<UInt>(supportStartIndices[0]); uInd != static_cast<UInt>(supportStartIndices[0]) + 3; ++uInd)
                indices.push_back(this->GetControlPointIndex(uInd, vInd));

        BSIPC_ASSERT(indices.size() == 9, "Support must have length 9");

        //std::sort(indices.begin(), indices.end());
        return indices;
    }

    /// Definitions of inline functions
    
    inline Float BSSurface::LengthSideX() const
    {
        BSIPC_ASSERT(!this->irregular, "LengthSideX only supports rectangle primitives");
        return (this->EvalAtRest(this->UDomainMax(), 0) - this->EvalAtRest(0, 0)).norm();
    }

    inline Float BSSurface::LengthSideY() const
    {
        BSIPC_ASSERT(!this->irregular, "LengthSideY only supports rectangle primitives");
        return (this->EvalAtRest(0, this->VDomainMax()) - this->EvalAtRest(0, 0)).norm();
    }

    inline Vec2 BSSurface::SingleDimSupport(UInt i, UInt bound) const
    {
        BSIPC_ASSERT(i >= 0 && i <= bound - 1, "Control point index i out of bound");

        Vec2 result = Vec2::Zero();
        if (i == 0)
            result = Vec2(0., 1.);
        else if (i == 1)
            result = Vec2(0., 2.);
        else if (i == bound - 2)
            result = Vec2(static_cast<Float>(i) - 2., static_cast<Float>(i));
        else if (i == bound - 1)
            result = Vec2(static_cast<Float>(i) - 2., static_cast<Float>(i) - 1.);
        else
            result = Vec2(static_cast<Float>(i) - 2., static_cast<Float>(i) + 1.);

        return result;
    }

    inline BSBasisFn BSSurface::SelectBasisFnByIndex(UInt i, const std::array<BSBasisFn, 5>& basis, UInt bound) const
    {
        BSBasisFn result = basis[2];
        if (i == 0)
            result = basis[0];
        else if (i == 1)
            result = basis[1];
        else if (i == bound - 2)
            result = basis[3];
        else if (i == bound - 1)
            result = basis[4];

        return result;
    }

    inline Bool BSSurface::CoordInNodeOpenSupport(UInt i, UInt j, Float u, Float v) const
    {
        BSIPC_ASSERT(i >= 0 && i <= this->resolutionX - 1, "Control Point i out of bound");
        BSIPC_ASSERT(j >= 0 && j <= this->resolutionY - 1, "Control Point j out of bound");
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        Vec2 uRange = this->SingleDimSupport(i, this->resolutionX);
        Vec2 vRange = this->SingleDimSupport(j, this->resolutionY);
       
        return InOpenInterval(u, uRange[0], uRange[1]) && InOpenInterval(v, vRange[0], vRange[1]);
    }

    inline Bool BSSurface::CoordInNodeClosedSupport(UInt i, UInt j, Float u, Float v) const
    {
        BSIPC_ASSERT(i >= 0 && i <= this->resolutionX - 1, "Control Point i out of bound");
        BSIPC_ASSERT(j >= 0 && j <= this->resolutionY - 1, "Control Point j out of bound");
        BSIPC_ASSERT(u >= 0 && u <= this->UDomainMax(), "U value out of bound");
        BSIPC_ASSERT(v >= 0 && v <= this->VDomainMax(), "V value out of bound");

        Vec2 uRange = this->SingleDimSupport(i, this->resolutionX);
        Vec2 vRange = this->SingleDimSupport(j, this->resolutionY);

        return InClosedInterval(u, uRange[0], uRange[1]) && InClosedInterval(v, vRange[0], vRange[1]);
    }
}
