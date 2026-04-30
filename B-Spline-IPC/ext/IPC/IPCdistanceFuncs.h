#pragma once
#include "Eigen/Eigen"
#include "mesh.h"
#include "collisionUtil.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

namespace IPC
{
    class BHessian {
    public:
        vector<int> D1Index;//pIndex, DpeIndex, DptIndex;
        vector<Vector3i> D3Index;
        vector<Vector4i> D4Index;
        vector<Vector2i> D2Index;//, DempIndex;
        vector<Matrix<double, 12, 12>> H12x12;//, HDee, HDpt;
        vector<Matrix<double, 3, 3>> H3x3;//, HDpm, HDpm, HDpm;
        vector<Matrix<double, 6, 6>> H6x6;//, HDeme;
        vector<Matrix<double, 9, 9>> H9x9;

        std::vector<Eigen::Triplet<double>> toTriplets(const vector<int>& Btype) const;
        void clear();
    };

    class IglUtils {
    public:
        template <typename Scalar, int size>
        static void makePD(Eigen::Matrix<Scalar, size, size>& symMtr)
        {
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, size, size>> eigenSolver(symMtr);
            if (eigenSolver.eigenvalues()[0] >= 0.0) {
                return;
            }
            Eigen::DiagonalMatrix<Scalar, size> D(eigenSolver.eigenvalues());
            int rows = ((size == Eigen::Dynamic) ? symMtr.rows() : size);
            for (int i = 0; i < rows; i++) {
                if (D.diagonal()[i] < 0.0) {
                    D.diagonal()[i] = 0.0;
                }
                else {
                    break;
                }
            }
            symMtr = eigenSolver.eigenvectors() * D * eigenSolver.eigenvectors().transpose();
        }

        static bool segTriIntersect(const Eigen::RowVector3d& ve0, const Eigen::RowVector3d& ve1,
            const Eigen::RowVector3d& vt0, const Eigen::RowVector3d& vt1, const Eigen::RowVector3d& vt2)
        {
            Eigen::Matrix3d coefMtr;
            coefMtr.col(0) = vt1 - vt0;
            coefMtr.col(1) = vt2 - vt0;
            coefMtr.col(2) = ve0 - ve1;

            Eigen::RowVector3d n = coefMtr.col(0).cross(coefMtr.col(1));
            if (n.dot(ve0 - vt0) * n.dot(ve1 - vt0) > 0.0) {
                return false; // edge is on one side of the plane that triangle is in
            }

            if (coefMtr.determinant() == 0.0) {
                return false; // coplanar, we can detect it by d(EE)=0 or d(PT)=0
            }

            Eigen::Vector3d uvt = coefMtr.fullPivLu().solve((ve0 - vt0).transpose());
            if (uvt[0] >= 0.0 && uvt[1] >= 0.0 && uvt[0] + uvt[1] <= 1.0 && uvt[2] >= 0.0 && uvt[2] <= 1.0) {
                return true;
            }
            else {
                return false;
            }
        }
    };

    int dType_PT(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        const Eigen::Vector3d& v3);

    int dType_EE(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        const Eigen::Vector3d& v3);

    void d_PP(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        double& d);

    void d_PE(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        double& d);

    void d_PT(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        const Eigen::Vector3d& v3,
        double& d);

    void computeEECrossSqNorm(
        const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        const Eigen::Vector3d& v3,
        double& result);

    void compute_e(
        const Eigen::RowVector3d& v0,
        const Eigen::RowVector3d& v1,
        const Eigen::RowVector3d& v2,
        const Eigen::RowVector3d& v3,
        double eps_x, double& e);

    void compute_eps_x(const mesh3D& mesh,
        int eI0, int eI1, int eJ0, int eJ1, double& eps_x);

    void d_EE(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        const Eigen::Vector3d& v3,
        double& d);
    void compute_b(double d, double dHat, double& b);
    void compute_g_b(double d, double dHat, double& g);
    void compute_H_b(double d, double dHat, double& H);

    int compute_g_dpt(const mesh3D& mesh,
        const std::vector<MMCVID>& activeSet,
        const Eigen::VectorXd& input,
        vector<Vector3d>& output_incremental,
        int offset,
        double coef);
    void compute_g_dpt_new(const mesh3D& mesh,
        const std::vector<MMCVID>& activeSet,
        const Eigen::VectorXd& input,
        vector<Vector3d>& output_incremental,
        int offset,
        double coef, double d_hat);
    void compute_g_dee(const mesh3D& mesh,
        //const std::vector<MMCVID>& paraEEMMCVIDSet,
        //const std::vector<std::pair<int, int>>& paraEEeIeJSet,
        vector<Vector3d>& grad_inc, double dHat, double coef);

    void Evaluate_SelfPTConstraintVals(const mesh3D& mesh, Eigen::VectorXd& vals, const int& offset);
    void Evaluate_SelfEEConstraintVals(const mesh3D& mesh, Eigen::VectorXd& vals, const int& offset);
    void Evaluate_GroundConstraintVals(const Ground& gd, const mesh3D& mesh, Eigen::VectorXd& vals, const int& offset);

    void compute_H_dpt(const mesh3D& mesh,
        BHessian& BH,
        double dHat, double coef);
    void compute_H_dpt_new(const mesh3D& mesh,
        BHessian& BH,
        double dHat, double coef, double d_hat);
    void compute_H_dee(const mesh3D& mesh,
        BHessian& BH,
        double dHat, double coef);

    void Self_largestFeasibleStepSize(const mesh3D& mesh,
        const SpatialHash& sh,
        const vector<Eigen::Vector3d>& searchDir,
        double slackness,
        std::vector<std::pair<int, int>>& candidates,
        double& stepSize);
    void Self_largestFeasibleStepSize_CCD(const mesh3D& mesh,
        SpatialHash& sh,
        const vector<Eigen::Vector3d>& searchDir,
        double slackness,
        double& stepSize);

    double SelfConstraintVal(const mesh3D& mesh, const MMCVID& active);

}
