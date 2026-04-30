#pragma once
#include "bsfwd.h"

#include "Geometry/BSpline/QuadPoint.h"
#include "Geometry/BSpline/BSSurface.h"

#include "Solver/SolverConfig.h"

namespace BSIPC
{
    namespace Energy
    {
        /* ***************************** Stretching ***************************** */

        // Uses FBW Energy

        /**** Stretching-Shearing Energy Conventions ****
         * Computing gradient and hessian of ss energy uses chain rule:
         * dPsi/dP = dF/dP * dPsi/dF: 3n-by-1 = 3n-by-6 * 6-by-1
         * d^2 Psi/dP^2 = (dF/dP)^t (d^2 Psi/dF^2) * (dF/dP): 3n-by-3n = 3n-by-6 * 6-by-6 * 6-by-3n
         *
         * Matrix flattening is done in {column-major} order (to coordinate with Hessian formula in FBW):
         * F = | F11 F12 | flattens to vec(F) = | F11 F21 F31 F12 F22 F32 |
         *     | F21 F22 |
         *     | F31 F32 |
         */
        
        Float SsValAt(const QuadPoint& quadPoint, const SolverConfig& config, const BSSurface* const surface);
        Float SsValAt(Mat<3, 2> F, const SolverConfig& config, const BSSurface* const surface);
        
        // Should return a 27-by-1 column vector (support is 9 nodes, each with 3 entries)
        DMat SsGradBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface);

        // Should return a 27-by-27 matrix (support is 9 nodes, each with 3 entries)
        Mat<27, 27> SsHessBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface);

        Float StVKSsValAt(const QuadPoint& quadPoint, const SolverConfig& config, const BSSurface* const surface);
        Float StVKSsValAt(Mat<3, 2> F, const SolverConfig& config, const BSSurface* const surface);
        Vec<27> StVKSsGradBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface);
        Mat<27, 27> StVKSsHessBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface);

#if defined BSIPC_NUMERIC_TEST 
        DMat NumericalSsGradBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, BSIPC_OUT BSSurface* const target);

        DMat NumericalSsHessBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, BSIPC_OUT BSSurface* const target);
#endif

        /* ***************************** Helper PDs ***************************** */

        // dPsi/dF, 6-by-1 matrix
        // @param PD_Param_Coord: dX/du and dX/dv evaluated at (u, v), where F is evaluated at
        DMat PD_SsEnergy_DeformGrad(const Mat<3, 2> F, const SolverConfig& config);

        // d^2Psi/dF^2, 6-by-6 matrix
        DMat PD2_SsEnergy_DeformGrad(const Mat<3, 2> F, const SolverConfig& config);

        // dF/dP_ij restricted to nonzero entries, 27-by-6 matrix
        DMat PD_DeformGrad_CtrlPt_Block(Float u, Float v, const std::vector<UInt> indices, const BSSurface* const surface);

        Vec<6> PD_StVKSsEnergy_DeformGrad(const Mat<3, 2> F, const SolverConfig& config);
        Mat<6, 6> PD2_StVKSsEnergy_DeformGrad(const Mat<3, 2> F, const SolverConfig& config);
    }
}