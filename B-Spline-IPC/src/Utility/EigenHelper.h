#pragma once

#include "bsfwd.h"

#include <iomanip>
#include <iostream>
#include <sstream>

#if defined BSIPC_USE_MKL
#include <Eigen/PardisoSupport>

#ifndef EIGEN_USE_MKL_ALL
#error "MKL is not enabled! EIGEN_USE_MKL_ALL not defined"
#endif
#endif

#include "Utility/Profiler.h"

// Simple helpers for Eigen-associated types. All functions should be inlined to avoid multiple definitions.

namespace BSIPC
{
    // Converting Eigen structures to String
    template <typename Scalar, Int Row, Int Col, Int Op, Int MaxRow, Int MaxCol>
    inline String ToStr(Eigen::Matrix<Scalar, Row, Col, Op, MaxRow, MaxCol> m)
    {
        std::stringstream ss;
        ss << "\n" << m;
        return ss.str();
    }

    inline String ToStr(SpMat m)
    {
        std::stringstream ss;
        ss << "\n" << m;
        return ss.str();
    }

    //// Convert Sparse Matrix from raw entry data
    //inline SpMat MakeSparseMat(UInt rowCnt, UInt colCnt, const SpMatData& entries)
    //{
    //    BSIPC_ASSERT(rowCnt * colCnt >= entries.size(), "Sparse matrix data has size larger than the specified matrix dimension");

    //    SpMat mat(rowCnt, colCnt);
    //    mat.setFromTriplets(entries.begin(), entries.end());
    //    return mat;
    //}

    // Convert Dynamic Matrix to Fixed Matrix
    template <Int X, Int Y>
    inline Mat<X, Y> FixateDynamicMatrix(DMat mat)
    {
        Int colNum = mat.cols();
        Int rowNum = mat.rows();
        BSIPC_ASSERT(rowNum <= X && colNum <= Y, "Resizing matrix with fewer dimensions will lost information");

        Mat<X, Y> fixedMat = mat.resize(X, Y);
        return fixedMat;
    }

    inline DMat SparseToDense(SpMat mat) { return DMat(mat); }

    inline SpMat DenseToSparse(DMat mat) { return DMat(mat).sparseView(); }

    template <Int X, Int Y>
    inline Mat<X, Y> SPDProjection(Mat<X, Y> mat)
    {
        Eigen::EigenSolver<DMat> es(mat);
        
        Eigen::VectorXd eigenvalues = es.eigenvalues().real();

        for (UInt i = 0; i != eigenvalues.size(); ++i)
            if (eigenvalues[i] < 0)
                eigenvalues[i] = 0;

        Eigen::MatrixXd V = es.eigenvectors().real();
        Eigen::MatrixXd Vt = V.transpose();
        Eigen::MatrixXd lambda = eigenvalues.asDiagonal();

        return V * lambda * Vt;
    }   

    inline DMat SolveLinearSystemEigenCG(const SpMat& A, const DMat& b)
    {
        BSIPC_ASSERT(A.rows() == b.rows(), "Matrix A and b have inconsistent dimensions");

        //Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<Float>> cg;
        Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<Float>> cg;

        BSIPC_INFO("A.size: {}x{}, B.size: {}x1", A.rows(), A.cols(), b.rows());
        BSIPC_INFO("nonzero entries in Hessian: {}", A.nonZeros());

        cg.setTolerance(1e-8);

        //BSIPC_PROFILE_SCOPE(
            //"AnalyzeCG",
            cg.analyzePattern(A);
        //)

        if (cg.info() != Eigen::Success)
            BSIPC_ERROR("AnalyzePattern failed");

        //BSIPC_PROFILE_SCOPE(
            //"factorizeCG",
            cg.factorize(A);
        //)

        if (cg.info() != Eigen::Success)
            BSIPC_ERROR("Factorize failed");

        //BSIPC_PROFILE_SCOPE(
            //"SolveCG",
            Eigen::VectorXd x = cg.solve(b);
        //)

        if (cg.info() != Eigen::Success)
            BSIPC_ERROR("Solve failed");

        //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > chol(A);
        //if (chol.info() != Eigen::Success)
        //    BSIPC_WARN("Cholesky decomposition failed");
        //Eigen::VectorXd x = chol.solve(b);
        //if (chol.info() != Eigen::Success)
        //    BSIPC_WARN("SimplicialLDLT solving failed");

        return x;
    }

    template<typename Scalar>
    bool ValidateSparseMatrix(const Eigen::SparseMatrix<Scalar>& mat, bool verbose = true) {
        bool is_valid = true;

        // Check if compressed
        if (!mat.isCompressed()) {
            if (verbose) std::cout << "ERROR: Matrix is not compressed. Call makeCompressed() first.\n";
            return false;
        }

        int rows = mat.rows();
        int cols = mat.cols();
        int nnz = mat.nonZeros();

        const int* inner_indices = mat.innerIndexPtr();
        const int* outer_starts = mat.outerIndexPtr();

        bool is_col_major = mat.IsRowMajor == 0;
        int outer_dim = is_col_major ? cols : rows;
        int inner_dim = is_col_major ? rows : cols;

        // Validate outer index pointer array
        if (outer_starts[0] != 0) {
            if (verbose) std::cout << "ERROR: outer_starts[0] = " << outer_starts[0] << ", expected 0\n";
            is_valid = false;
        }

        if (outer_starts[outer_dim] != nnz) {
            if (verbose) std::cout << "ERROR: outer_starts[" << outer_dim << "] = " << outer_starts[outer_dim]
                << ", expected " << nnz << "\n";
            is_valid = false;
        }

        // Check outer indices are non-decreasing
        for (int i = 0; i < outer_dim; ++i) {
            if (outer_starts[i] > outer_starts[i + 1]) {
                if (verbose) std::cout << "ERROR: outer_starts[" << i << "] = " << outer_starts[i]
                    << " > outer_starts[" << i + 1 << "] = " << outer_starts[i + 1] << "\n";
                is_valid = false;
            }

            // Validate inner indices for this outer index
            int start = outer_starts[i];
            int end = outer_starts[i + 1];

            for (int j = start; j < end; ++j) {
                // Check bounds
                if (inner_indices[j] < 0 || inner_indices[j] >= inner_dim) {
                    if (verbose) std::cout << "ERROR: inner_indices[" << j << "] = " << inner_indices[j]
                        << " out of range [0, " << inner_dim << ")\n";
                    is_valid = false;
                }

                // Check sorted order within each outer index
                if (j > start && inner_indices[j] <= inner_indices[j - 1]) {
                    if (verbose) std::cout << "ERROR: inner_indices not sorted at outer index " << i
                        << ": inner_indices[" << j - 1 << "] = " << inner_indices[j - 1]
                        << " >= inner_indices[" << j << "] = " << inner_indices[j] << "\n";
                    is_valid = false;
                }
            }
        }

        if (verbose && is_valid) {
            std::cout << "  Sparse matrix indices are VALID\n";
            std::cout << "  Format: " << (is_col_major ? "Column-major (CSC)" : "Row-major (CSR)") << "\n";
            std::cout << "  Dimensions: " << rows << " x " << cols << "\n";
            std::cout << "  Non-zeros: " << nnz << "\n";
        }

        return is_valid;
    }

    inline Bool SolveLinearSystemLDLT(const SpMat& A, const DMat& b, BSIPC_OUT DMat& x)
    {
        BSIPC_ASSERT(A.rows() == b.rows(), "Matrix A and b have inconsistent dimensions");

        Bool success = true;

        BSIPC_PROFILE_START("SolveLinearSystemLDLT");
#if defined BSIPC_USE_MKL
        Eigen::PardisoLU<SpMat> ldlt;
        //ldlt.pardisoParameterArray()[10] = 1;
        //ldlt.pardisoParameterArray()[12] = 1;
        //ldlt.pardisoParameterArray()[20] = 1;
        BSIPC_INFO("Using MKL PardisoLDLT solver");
#else
        Eigen::SimplicialLDLT<SpMat> ldlt;
#endif
        ldlt.compute(A);
        if (ldlt.info() != Eigen::Success)
            BSIPC_ERROR("LDLT factorization failed");

        x = ldlt.solve(b);
        BSIPC_ASSERT(x.rows() == b.rows(), "Solve result should have same shape as gradient");

        if (ldlt.info() != Eigen::Success)
            BSIPC_ERROR("LDLT solve failed");
        BSIPC_PROFILE_END("SolveLinearSystemLDLT");

        DMat diff = A * x - b;
        BSIPC_PEEK(diff.lpNorm<Eigen::Infinity>());
        BSIPC_PEEK(b.lpNorm<Eigen::Infinity>());
        if (diff.lpNorm<Eigen::Infinity>() > 0.0001 * b.lpNorm<Eigen::Infinity>())
        {
            BSIPC_WARN("LDLT solution verification failed with residual norm {} compared with RHS norm {}", diff.lpNorm<Eigen::Infinity>(), b.lpNorm<Eigen::Infinity>());
            Eigen::SimplicialLDLT<SpMat> ldlt;
            ldlt.compute(A);
            if (ldlt.info() != Eigen::Success)
                BSIPC_ERROR("SimpliclalDLT factorization failed on re-compute");
            x = ldlt.solve(b);
            if (ldlt.info() != Eigen::Success)
                BSIPC_ERROR("SimpliclalDLT solve failed on re-compute");
            diff = A * x - b;

            if (diff.lpNorm<Eigen::Infinity>() > 0.001 * b.lpNorm<Eigen::Infinity>())
            {
                BSIPC_ERROR(
                    "LDLT solution verification failed again after re-compute with residual norm {} compared with RHS norm {}", 
                    diff.lpNorm<Eigen::Infinity>(), b.lpNorm<Eigen::Infinity>()
                );
                success = false;
            }
        }

        return success;
    }

    inline DMat SPDProjection(DMat mat)
    {
        Eigen::EigenSolver<DMat> es(mat);
        
        Eigen::VectorXd eigenvalues = es.eigenvalues().real();

        for (UInt i = 0; i != eigenvalues.size(); ++i)
            if (eigenvalues[i] < 0)
                eigenvalues[i] = 0;

        Eigen::MatrixXd V = es.eigenvectors().real();
        Eigen::MatrixXd Vt = V.transpose();
        Eigen::MatrixXd lambda = eigenvalues.asDiagonal();

        return V * lambda * Vt;
    }

    // returns A * B
    inline void SpMatProduct(const SpMat& A, const SpMat& B, BSIPC_OUT SpMat& result)
    {
#if defined BSIPC_USE_MKL
        sparse_matrix_t mkl_A, mkl_B, mkl_C;

        mkl_sparse_d_create_csc(&mkl_A, SPARSE_INDEX_BASE_ZERO,
            A.rows(), A.cols(),
            const_cast<Int*>(A.outerIndexPtr()),
            const_cast<Int*>(A.outerIndexPtr() + 1),
            const_cast<Int*>(A.innerIndexPtr()),
            const_cast<Float*>(A.valuePtr()));

        mkl_sparse_d_create_csc(&mkl_B, SPARSE_INDEX_BASE_ZERO,
            B.rows(), B.cols(),
            const_cast<Int*>(B.outerIndexPtr()),
            const_cast<Int*>(B.outerIndexPtr() + 1),
            const_cast<Int*>(B.innerIndexPtr()),
            const_cast<Float*>(B.valuePtr()));

        // Multiply: C = A * B
        sparse_status_t status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, mkl_A, mkl_B, &mkl_C);
        if (status != SPARSE_STATUS_SUCCESS)
            BSIPC_ERROR("MKL sparse matrix multiplication failed with status code {}", static_cast<UInt>(status));

        sparse_index_base_t indexing;
        MKL_INT rows, cols;
        MKL_INT* cols_start, * cols_end, * row_indx;
        Float* values;

        mkl_sparse_d_export_csc(mkl_C, &indexing, &rows, &cols,
            &cols_start, &cols_end, &row_indx, &values);

        MKL_INT nnz = cols_end[cols - 1];

        result.resize(rows, cols);
        result.resizeNonZeros(nnz);

        std::memcpy(result.outerIndexPtr(), cols_start, (cols + 1) * sizeof(Int));
        std::memcpy(result.innerIndexPtr(), row_indx, nnz * sizeof(Int));
        std::memcpy(result.valuePtr(), values, nnz * sizeof(Float));

        result.makeCompressed();

        // Cleanup
        mkl_sparse_destroy(mkl_A);
        mkl_sparse_destroy(mkl_B);
        mkl_sparse_destroy(mkl_C);
#else
        result = A * B;
#endif
    }

    // computes proj^T * perm * hess * perm^T * proj
    inline void ProjectKKTHessian(const SpMat& proj, const SpMat& perm, const SpMat& hess, BSIPC_OUT SpMat& result)
    {
#if defined BSIPC_USE_MKL
        sparse_matrix_t mkl_proj, mkl_perm, mkl_projPerm, mkl_projPerm_T, mkl_hess, mkl_temp, mkl_result;

        // Create MKL sparse matrices
        mkl_sparse_d_create_csc(&mkl_proj, SPARSE_INDEX_BASE_ZERO,
            proj.rows(), proj.cols(),
            const_cast<Int*>(proj.outerIndexPtr()),
            const_cast<Int*>(proj.outerIndexPtr() + 1),
            const_cast<Int*>(proj.innerIndexPtr()),
            const_cast<Float*>(proj.valuePtr()));

        mkl_sparse_d_create_csc(&mkl_perm, SPARSE_INDEX_BASE_ZERO,
            perm.rows(), perm.cols(),
            const_cast<Int*>(perm.outerIndexPtr()),
            const_cast<Int*>(perm.outerIndexPtr() + 1),
            const_cast<Int*>(perm.innerIndexPtr()),
            const_cast<Float*>(perm.valuePtr()));

        mkl_sparse_d_create_csc(&mkl_hess, SPARSE_INDEX_BASE_ZERO,
            hess.rows(), hess.cols(),
            const_cast<Int*>(hess.outerIndexPtr()),
            const_cast<Int*>(hess.outerIndexPtr() + 1),
            const_cast<Int*>(hess.innerIndexPtr()),
            const_cast<Float*>(hess.valuePtr()));

        BSIPC_PROFILE_START("compute spmm");
        // Compute proj^T * perm
        sparse_status_t status = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, mkl_proj, mkl_perm, &mkl_projPerm);
        if (status != SPARSE_STATUS_SUCCESS)
            BSIPC_ERROR("MKL sparse matrix multiplication (proj^T * perm) failed with status code {}", static_cast<UInt>(status));

        // Compute projPerm^T
        status = mkl_sparse_convert_csc(mkl_projPerm, SPARSE_OPERATION_TRANSPOSE, &mkl_projPerm_T);
        if (status != SPARSE_STATUS_SUCCESS)
            BSIPC_ERROR("MKL sparse matrix transpose (projPerm^T) failed with status code {}", static_cast<UInt>(status));

        // Compute temp = projPerm * hess
        status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, mkl_projPerm, mkl_hess, &mkl_temp);
        if (status != SPARSE_STATUS_SUCCESS)
            BSIPC_ERROR("MKL sparse matrix multiplication (projPerm * hess) failed with status code {}", static_cast<UInt>(status));

        // Compute result = temp * projPerm^T
        status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, mkl_temp, mkl_projPerm_T, &mkl_result);
        if (status != SPARSE_STATUS_SUCCESS)
            BSIPC_ERROR("MKL sparse matrix multiplication (temp * projPerm^T) failed with status code {}", static_cast<UInt>(status));

        mkl_sparse_order(mkl_result);
        BSIPC_PROFILE_END("compute spmm");

        // Convert back
        sparse_index_base_t indexing;
        MKL_INT rows, cols;
        MKL_INT* cols_start, * cols_end, * row_indx;
        Float* values;

        mkl_sparse_d_export_csc(mkl_result, &indexing, &rows, &cols,
            &cols_start, &cols_end, &row_indx, &values);

        MKL_INT nnz = cols_end[cols - 1];

        result.resize(rows, cols);
        result.resizeNonZeros(nnz);

        std::memcpy(result.outerIndexPtr(), cols_start, (cols + 1) * sizeof(Int));
        std::memcpy(result.innerIndexPtr(), row_indx, nnz * sizeof(Int));
        std::memcpy(result.valuePtr(), values, nnz * sizeof(Float));

        result.makeCompressed();

        // Cleanup
        mkl_sparse_destroy(mkl_proj);
        mkl_sparse_destroy(mkl_perm);
        mkl_sparse_destroy(mkl_projPerm);
        mkl_sparse_destroy(mkl_projPerm_T);
        mkl_sparse_destroy(mkl_hess);
        mkl_sparse_destroy(mkl_temp);
        mkl_sparse_destroy(mkl_result);
#else
        SpMat projPerm = proj.transpose() * perm;
        result = projPerm * hess * projPerm.transpose();
#endif
    }

    template <typename SpMat>
    void PrintSparseMatrixStructure(const SpMat& M, const std::string& name = "SparseMatrix")
    {
        using Scalar = typename SpMat::Scalar;
        using Index = typename SpMat::StorageIndex;

        std::cout << "==== " << name << " ====\n";
        std::cout << "Dimensions: " << M.rows() << " x " << M.cols()
            << ", nnz = " << M.nonZeros() << "\n";
        std::cout << "isCompressed = " << M.isCompressed() << "\n";

        // --- Outer indices ---
        std::cout << "outerIndexPtr: { ";
        if (M.outerIndexPtr() != nullptr) {
            for (int i = 0; i < M.outerSize() + 1; ++i)
                std::cout << M.outerIndexPtr()[i] << (i < M.outerSize() ? ", " : "");
        }
        else std::cout << "(null)";
        std::cout << " }\n";

        // --- Inner indices ---
        std::cout << "innerIndexPtr: { ";
        if (M.innerIndexPtr() != nullptr) {
            for (int i = 0; i < M.nonZeros(); ++i)
                std::cout << M.innerIndexPtr()[i] << (i < M.nonZeros() - 1 ? ", " : "");
        }
        else std::cout << "(null)";
        std::cout << " }\n";

        // --- Values ---
        std::cout << "valuePtr: { ";
        if (M.valuePtr() != nullptr) {
            std::cout << std::fixed << std::setprecision(6);
            for (int i = 0; i < M.nonZeros(); ++i)
                std::cout << M.valuePtr()[i] << (i < M.nonZeros() - 1 ? ", " : "");
        }
        else std::cout << "(null)";
        std::cout << " }\n";

        std::cout << "===========================\n";
    }
}
