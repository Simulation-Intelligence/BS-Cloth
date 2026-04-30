#pragma once
#include "bsfwd.h"

#include "suitesparse/cholmod/cholmod.h"

namespace BSIPC
{
    class CholmodSolver
    {
    private:
        cholmod_common c;
        cholmod_sparse A_cholmod;
        cholmod_dense* b_cholmod;
        cholmod_factor* L;
        cholmod_dense* x_cholmod;

    public:
        CholmodSolver();
        ~CholmodSolver();

        DMat SolveLinearSystem(SpMat& A, const Eigen::VectorXd& b);
    };
}

