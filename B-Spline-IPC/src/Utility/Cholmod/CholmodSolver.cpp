#include "bspch.h"
#include "CholmodSolver.h"

#include "Utility/Profiler.h"

namespace BSIPC
{
    CholmodSolver::CholmodSolver()
    {
        cholmod_start(&c);

        A_cholmod.z = nullptr;
        A_cholmod.stype = 1; // symmetric upper part
        A_cholmod.itype = CHOLMOD_INT;
        A_cholmod.xtype = CHOLMOD_REAL;
        A_cholmod.dtype = CHOLMOD_DOUBLE;
        A_cholmod.sorted = 1;
        A_cholmod.packed = 1;

        b_cholmod = nullptr;
        x_cholmod = nullptr;
        L = nullptr;
    }

    CholmodSolver::~CholmodSolver()
    {
        // Clean up
        cholmod_free_dense(&x_cholmod, &c);
        cholmod_free_dense(&b_cholmod, &c);
        cholmod_free_factor(&L, &c);
        cholmod_finish(&c);
    }

    DMat CholmodSolver::SolveLinearSystem(SpMat& A, const Eigen::VectorXd& b)
    {
        // Convert Eigen::SparseMatrix to cholmod_sparse
        A.makeCompressed();

        A_cholmod.nrow = A.rows();
        A_cholmod.ncol = A.cols();
        A_cholmod.nzmax = A.nonZeros();
        A_cholmod.p = const_cast<int*>(A.outerIndexPtr());
        A_cholmod.i = const_cast<int*>(A.innerIndexPtr());
        A_cholmod.x = const_cast<double*>(A.valuePtr());

        // Convert RHS vector b to cholmod_dense
        if (b_cholmod == nullptr)
            b_cholmod = cholmod_allocate_dense(b.size(), 1, b.size(), CHOLMOD_REAL, &c);

        std::memcpy(b_cholmod->x, b.data(), b.size() * sizeof(double));

        BSIPC_PROFILE_START("analyze");

        // Factorize
        L = cholmod_analyze(&A_cholmod, &c);
        cholmod_factorize(&A_cholmod, L, &c);

        BSIPC_PROFILE_END("analyze");
        BSIPC_PROFILE_START("solve");

        // Solve
        x_cholmod = cholmod_solve(CHOLMOD_A, L, b_cholmod, &c);

        BSIPC_PROFILE_END("solve");

        // Copy result back to Eigen
        double* x_vals = static_cast<double*>(x_cholmod->x);
        DMat result = Eigen::Map<Eigen::VectorXd>(x_vals, b.size());

        return result;
    }
}
