#pragma once
#include "bsfwd.h"

#include "Geometry/BSpline/BSSurface.h"
#include "Geometry/BSpline/QuadPoint.h"

#include "Solver/BSTargetInfo.h"
#include "Solver/SolverConfig.h"

namespace BSIPC
{
	namespace Energy
	{
        /* ***************************** Bending ***************************** */

        // Uses Quadratic Bending Model (computed directly from mean curvature)

        // Bending energy density at (u, v)
        Float BdValAt(const QuadPoint& quadPoint, const SolverConfig& config, const BSSurface* const surface);

        // Should return a 27-by-1 column vector (support is 9 nodes, each with 3 entries)
        DMat BdGradBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface);

        // Should return a 27-by-27 matrix (support is 9 nodes, each with 3 entries)
        Mat<27, 27> BdHessBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface);

#if defined BSIPC_NUMERIC_TEST
        DMat NumericalBdGradBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, BSIPC_OUT BSSurface* const target);

        DMat NumericalBdHessBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, BSIPC_OUT BSSurface* const target);
        DMat NumericalBdHess(const SolverConfig& config, const std::vector<BSTargetInfo>& infos);
#endif
	}
}