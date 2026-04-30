#include "bspch.h"

#include "Membrane.h"

#include "Utility/EigenHelper.h"

namespace BSIPC
{
    namespace Energy
    {
        // TODO cache PD_Param_Coord_RestSpace
        Float SsValAt(const QuadPoint& quadPoint, const SolverConfig& config, const BSSurface* const surface)
        {
            Mat<3, 2> F = quadPoint.deformationGradient;

            Float result = SsValAt(F, config, surface);

            return result;
        }

        Float SsValAt(Mat<3, 2> F, const SolverConfig& config, const BSSurface* const surface)
        {
            Vec2 coordX1Dir = Vec2(1, 0);
            Vec2 coordX2Dir = Vec2(0, 1);

            Float I6 = coordX1Dir.transpose() * F.transpose() * F * coordX2Dir;
            Float shearEnergy = I6 * I6;

            Float I5u = (F * coordX1Dir).norm();
            Float I5v = (F * coordX2Dir).norm();

            Float uCoeff = 1;
            Float vCoeff = 1;

            if (I5u <= 1)
                uCoeff = 0;
            if (I5v <= 1)
                vCoeff = 0;

            Float stretchEnergy =
                pow(I5u - 1, 2) + uCoeff * config.strainRate * pow(I5u - 1, 3) +
                pow(I5v - 1, 2) + vCoeff * config.strainRate * pow(I5v - 1, 3);

            Float result = (config.stretchStiffness * stretchEnergy + config.shearStiffness * shearEnergy);

            return config.thickness * result;
        }

        DMat SsGradBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface)
        {
            DMat mat = DMat::Zero(27, 1);       // support contains 9 nodes, each with 3 coordinate entries

            const Mat<3, 2>& F = quadPoint.deformationGradient;

            //this->AppendToLog(ToStr(F));

            DMat pd_SsEnergy_DeformGrad = PD_SsEnergy_DeformGrad(F, config);

            DMat pd_DeformGrad_CtrlPt = PD_DeformGrad_CtrlPt_Block(quadPoint.U(), quadPoint.V(), indices, surface);

            mat = pd_DeformGrad_CtrlPt * pd_SsEnergy_DeformGrad;

            return config.thickness * mat;
        }

        Mat<27, 27> SsHessBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface)
        {
            Mat<27, 27> mat = Mat<27, 27>::Zero();           // support contains 9 nodes, each with 3 coordinate entries

            const Mat<3, 2>& F = quadPoint.deformationGradient;

            DMat pd2_SsEnergy_DeformGrad = PD2_SsEnergy_DeformGrad(F, config);
            DMat pd_DeformGrad_CtrlPt = PD_DeformGrad_CtrlPt_Block(quadPoint.U(), quadPoint.V(), indices, surface);

            // TODO test whether we can fill this in manually
            mat = pd_DeformGrad_CtrlPt * pd2_SsEnergy_DeformGrad * pd_DeformGrad_CtrlPt.transpose();
            Mat<27, 27> hess = config.thickness * mat;

            return hess;
        }

#if defined BSIPC_NUMERIC_TEST
        DMat NumericalSsGradBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, BSIPC_OUT BSSurface* const target)
        {
            DMat grad = DMat::Zero(3 * indices.size(), 1);
            const Float h = 1e-8;

            for (UInt i = 0; i != indices.size(); ++i)
            {
                for (UInt index = 0; index != 3; ++index)
                {
                    UInt nodeIndex = indices[i];
                    Vec3& perturbVertex = target->GetControlPoint(nodeIndex);

                    // Up, f(... x_i + h ...)
                    perturbVertex[index] += h;
                    target->UpdateQuadPointCache();
                    Float upVal = SsValAt(quadPoint, config, target);

                    // Down, f(... x_i - h ...)
                    perturbVertex[index] -= 2 * h;
                    target->UpdateQuadPointCache();
                    Float downVal = SsValAt(quadPoint, config, target);

                    perturbVertex[index] += h;

                    grad(3 * i + index, 0) = (upVal - downVal) / (2 * h);
                }
            }

            return grad;
        }

        DMat NumericalSsHessBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, BSIPC_OUT BSSurface* const target)
        {
            const Float h = 1e-8;
            const Float hSq = h * h;
            DMat hess = DMat::Zero(27, 27);

            for (UInt i = 0; i != indices.size(); ++i)
            {
                for (UInt j = 0; j != indices.size(); ++j)
                {
                    for (UInt indexI = 0; indexI != 3; ++indexI)
                    {
                        for (UInt indexJ = 0; indexJ != 3; ++indexJ)
                        {
                            UInt matIndexI = 3 * i + indexI;
                            UInt matIndexJ = 3 * j + indexJ;

                            Vec3& perturbVertexI = target->GetControlPoint(indices[i]);
                            Vec3& perturbVertexJ = target->GetControlPoint(indices[j]);

                            // Up, f(... x_i + h ... x_j + h ...)
                            perturbVertexI[indexI] += h;
                            perturbVertexJ[indexJ] += h;
                            target->UpdateQuadPointCache();
                            Float upSsVal = SsValAt(quadPoint, config, target);

                            // Left, f(... x_i - h ... x_j + h ...)
                            perturbVertexI[indexI] -= 2 * h;
                            target->UpdateQuadPointCache();
                            Float leftSsVal = SsValAt(quadPoint, config, target);

                            // Right, f(... x_i + h ... x_j - h ...)
                            perturbVertexI[indexI] += 2 * h;
                            perturbVertexJ[indexJ] -= 2 * h;
                            target->UpdateQuadPointCache();
                            Float rightSsVal = SsValAt(quadPoint, config, target);

                            // Down, f(... x_i - h ... x_j - h ...)
                            perturbVertexI[indexI] -= 2 * h;
                            target->UpdateQuadPointCache();
                            Float downSsVal = SsValAt(quadPoint, config, target);

                            Float hessEntry = (upSsVal - rightSsVal - leftSsVal + downSsVal) / (4 * hSq);
                            hess(matIndexI, matIndexJ) = hessEntry;
                        }
                    }
                }
            }

            return hess;
        }
#endif

        DMat PD_SsEnergy_DeformGrad(const Mat<3, 2> F, const SolverConfig& config)
        {
            Vec2 coordX1Dir = Vec2(1, 0);
            Vec2 coordX2Dir = Vec2(0, 1);

            Float I6 = coordX1Dir.transpose() * F.transpose() * F * coordX2Dir;
            Mat<3, 2> stretchPk1, shearPk1;

            shearPk1 = 2 * (I6 - coordX1Dir.transpose() * coordX2Dir) * (
                F * coordX1Dir * coordX2Dir.transpose() +
                F * coordX2Dir * coordX1Dir.transpose()
                );

            Float I5u = (F * coordX1Dir).transpose() * F * coordX1Dir;
            Float I5v = (F * coordX2Dir).transpose() * F * coordX2Dir;

            Float uCoeff = 1.0 - 1 / sqrt(I5u);
            Float vCoeff = 1.0 - 1 / sqrt(I5v);

            if (I5u > 1)
                uCoeff += 1.5 * config.strainRate * (sqrt(I5u) + 1 / sqrt(I5u) - 2);
            if (I5v > 1)
                vCoeff += 1.5 * config.strainRate * (sqrt(I5v) + 1 / sqrt(I5v) - 2);

            stretchPk1 = uCoeff * 2. * F * coordX1Dir * coordX1Dir.transpose() + vCoeff * 2. * F * coordX2Dir * coordX2Dir.transpose();
            Mat<3, 2> result = config.stretchStiffness * stretchPk1 + config.shearStiffness * shearPk1;

            DMat reshaped;
            reshaped.resize(6, 1);
            reshaped <<
                result(0, 0), result(1, 0), result(2, 0), result(0, 1), result(1, 1), result(2, 1);

#if defined BSIPC_NUMERIC_TEST
            //Float du = 1e-6;
            //DMat reshapedNumerical = DMat::Zero(6, 1);
            //for (UInt i = 0; i != 6; ++i)
            //{
            //    Mat<3, 2> perturbF = F;
            //    perturbF(i / 2, i % 2) += du;
            //    Float perturbEnergy = this->SsValAt(perturbF);
            //    reshapedNumerical(i, 0) = (perturbEnergy - this->SsValAt(F)) / du;
            //}

            //this->AppendToLog("===== PD SsEnergy DeformGrad Test =====\n\n");
            //this->AppendToLog("Analytic PD SsEnergy DeformGrad: \n" + ToStr(reshaped) + "\n");
            //this->AppendToLog("Numerical PD SsEnergy DeformGrad: \n" + ToStr(reshapedNumerical) + "\n");
            //this->AppendToLog("\n<<<<< End Test >>>>>\n");
#endif

            return reshaped;
        }

        DMat PD2_SsEnergy_DeformGrad(const Mat<3, 2> F, const SolverConfig& config)
        {
            Mat<6, 6> finalH = Mat<6, 6>::Zero();

            Vec2 coordX1Dir = Vec2(1, 0);
            Vec2 coordX2Dir = Vec2(0, 1);

            {
                Eigen::Matrix<double, 6, 6> H;
                H.setZero();
                double I5u = (F * coordX1Dir).transpose() * F * coordX1Dir;
                double I5v = (F * coordX2Dir).transpose() * F * coordX2Dir;

                double sqrtI5u = sqrt(I5u);
                double sqrtI5v = sqrt(I5v);

                if (sqrtI5u > 1)
                    H(0, 0) = H(1, 1) = H(2, 2) = 2 * (((sqrtI5u - 1) * (3 * sqrtI5u * config.strainRate - 3 * config.strainRate + 2)) / (2 * sqrtI5u));
                if (sqrtI5v > 1)
                    H(3, 3) = H(4, 4) = H(5, 5) = 2 * (((sqrtI5v - 1) * (3 * sqrtI5v * config.strainRate - 3 * config.strainRate + 2)) / (2 * sqrtI5v));
                auto fu = F.col(0).normalized();
                auto fv = F.col(1).normalized();

                double uCoeff = (sqrtI5u > 1.0) ? (3 * I5u * config.strainRate - 3 * config.strainRate + 2)
                    / (sqrt(I5u)) :
                    2.0;
                double vCoeff = (sqrtI5v > 1.0) ? (3 * I5v * config.strainRate - 3 * config.strainRate + 2)
                    / (sqrt(I5v)) :
                    2.0;


                H.block<3, 3>(0, 0) += uCoeff * (fu * fu.transpose());
                H.block<3, 3>(3, 3) += vCoeff * (fv * fv.transpose());

                finalH += config.stretchStiffness * H;
            }
            {
                Eigen::Matrix<double, 6, 6> H_shear;
                H_shear.setZero();
                Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
                H(3, 0) = H(4, 1) =
                    H(5, 2) = H(0, 3) =
                    H(1, 4) = H(2, 5) = 1.0;
                double I6 = coordX1Dir.transpose() * F.transpose() * F * coordX2Dir;
                double signI6 = (I6 >= 0) ? 1.0 : -1.0;
                auto g = F * (coordX1Dir * coordX2Dir.transpose() + coordX2Dir * coordX1Dir.transpose());
                Eigen::Matrix<double, 6, 1> vec_g = Eigen::Matrix<double, 6, 1>::Zero();

                vec_g.block(0, 0, 3, 1) = g.col(0);
                vec_g.block(3, 0, 3, 1) = g.col(1);
                double I2 = F.squaredNorm();
                double lambda0 = 0.5 * (I2 + sqrt(I2 * I2 + 12.0 * I6 * I6));
                Eigen::Matrix<double, 6, 1> q0 = (I6 * H * vec_g + lambda0 * vec_g).normalized();
                Eigen::Matrix<double, 6, 6> T = Eigen::Matrix<double, 6, 6>::Identity();
                T = 0.5 * (T + signI6 * H);
                auto Tq = T * q0;
                double normTq = Tq.squaredNorm();
                H_shear = fabs(I6) * (T - (Tq * Tq.transpose()) / normTq) + lambda0 * (q0 * q0.transpose());
                H_shear *= 2;
                finalH += config.shearStiffness * H_shear;
            }

            // The followings are non-SPD stretching local Hessian
            //{
            //    Eigen::Matrix<double, 6, 6> H;
            //    H.setZero();
            //    double I5u = (F * coordX1Dir).transpose() * F * coordX1Dir;
            //    double I5v = (F * coordX2Dir).transpose() * F * coordX2Dir;

            //    double sqrtI5u = sqrt(I5u);
            //    double sqrtI5v = sqrt(I5v);

            //    H(0, 0) = H(1, 1) = H(2, 2) = 2 * (((sqrtI5u - 1) * (3 * sqrtI5u * config.strainRate - 3 * config.strainRate + 2)) / (2 * sqrtI5u));
            //    H(3, 3) = H(4, 4) = H(5, 5) = 2 * (((sqrtI5v - 1) * (3 * sqrtI5v * config.strainRate - 3 * config.strainRate + 2)) / (2 * sqrtI5v));
            //    auto fu = F.col(0).normalized();
            //    auto fv = F.col(1).normalized();

            //    double uCoeff = (3 * I5u * config.strainRate - 3 * config.strainRate + 2)
            //        / (sqrt(I5u));
            //    double vCoeff = (3 * I5v * config.strainRate - 3 * config.strainRate + 2)
            //        / (sqrt(I5v));


            //    H.block<3, 3>(0, 0) += uCoeff * (fu * fu.transpose());
            //    H.block<3, 3>(3, 3) += vCoeff * (fv * fv.transpose());

            //    finalH += config.stretchStiffness * H;
            //}
            //{
            //    Eigen::Matrix<double, 6, 6> H_shear;
            //    H_shear.setZero();
            //    Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
            //    H(3, 0) = H(4, 1) =
            //        H(5, 2) = H(0, 3) =
            //        H(1, 4) = H(2, 5) = 1.0;
            //    double I6 = coordX1Dir.transpose() * F.transpose() * F * coordX2Dir;
            //    double signI6 = 1.;
            //    auto g = F * (coordX1Dir * coordX2Dir.transpose() + coordX2Dir * coordX1Dir.transpose());
            //    Eigen::Matrix<double, 6, 1> vec_g = Eigen::Matrix<double, 6, 1>::Zero();

            //    vec_g.block(0, 0, 3, 1) = g.col(0);
            //    vec_g.block(3, 0, 3, 1) = g.col(1);
            //    double I2 = F.squaredNorm();
            //    double lambda0 = 0.5 * (I2 + sqrt(I2 * I2 + 12.0 * I6 * I6));
            //    Eigen::Matrix<double, 6, 1> q0 = (I6 * H * vec_g + lambda0 * vec_g).normalized();
            //    Eigen::Matrix<double, 6, 6> T = Eigen::Matrix<double, 6, 6>::Identity();
            //    T = 0.5 * (T + signI6 * H);
            //    auto Tq = T * q0;
            //    double normTq = Tq.squaredNorm();
            //    H_shear = (I6) * (T - (Tq * Tq.transpose()) / normTq) + lambda0 * (q0 * q0.transpose());
            //    H_shear *= 2;
            //    finalH += config.shearStiffness * H_shear;
            //}

            return finalH;
        }

        DMat PD_DeformGrad_CtrlPt_Block(Float u, Float v, const std::vector<UInt> indices, const BSSurface* const target)
        {
            DMat mat = DMat::Zero(27, 6);

            for (UInt ptInd = 0; ptInd != indices.size(); ++ptInd)
            {
                IVec2 indPair = target->GetControlPointPosIndex(indices[ptInd]);
                UInt i = indPair[0], j = indPair[1];

                Vec2 pd = target->PD_DeformGrad_CtrlPt_At(i, j, u, v);

                for (UInt coordInd = 0; coordInd != 3; ++coordInd)
                {
                    mat(3 * ptInd + coordInd, coordInd) = pd[0];
                    mat(3 * ptInd + coordInd, coordInd + 3) = pd[1];
                }
            }

            return mat;
        }

        Float StVKSsValAt(const QuadPoint& quadPoint, const SolverConfig& config, const BSSurface* const surface)
        {
            Mat<3, 2> F = quadPoint.deformationGradient;

            Float result = StVKSsValAt(F, config, surface);

            return result;
        }

        // 1. Energy value
        Float StVKSsValAt(Mat<3, 2> F, const SolverConfig& config, const BSSurface* const surface)
        {
            Float lambda = config.lambda;
            Float mu = config.mu;
            Float thickness = config.thickness;

            Mat<2, 2> E = 0.5 * (F.transpose() * F - Mat<2, 2>::Identity());
            Float energy = 0.5 * lambda * E.trace() * E.trace() + mu * (E * E).trace();

            return thickness * energy;
        }

        // 2. First derivative (gradient) w.r.t. F
        Vec<6> PD_StVKSsEnergy_DeformGrad(const Mat<3, 2> F, const SolverConfig& config)
        {
            Vec<6> grad = Vec<6>::Zero();

            Float lambda = config.lambda;
            Float mu = config.mu;
            Float thickness = config.thickness;

            const Mat<2, 2> E = 0.5 * (F.transpose() * F - Mat<2, 2>::Identity());
            const auto traceE = E.trace();
            const Float psi = mu * E.squaredNorm() + 0.5 * lambda* traceE * traceE;

            const Mat<2, 2> dPsi_dE = 2.0 * mu * E + lambda * traceE * Mat<2, 2>::Identity();
            Eigen::Map<Mat<3, 2>> gPsi_dF(grad.data());
            gPsi_dF = F * dPsi_dE;

            grad *= thickness;

            return grad;
        }

        inline void svd_2x2(const Mat<2, 2>& A, Mat<2, 2>& U, Vec2& Sigma, Mat<2, 2>& V) {
            // Compute Su = A * A'
            Mat<2, 2> Su = A * A.transpose();

            // Find phi
            Float phi = 0.5 * atan2(Su(0, 1) + Su(1, 0), Su(0, 0) - Su(1, 1));
            Float Cphi = cos(phi);
            Float Sphi = sin(phi);
            U << Cphi, -Sphi, Sphi, Cphi;

            // Compute Sw = A' * A
            Mat<2, 2> Sw = A.transpose() * A;

            // Find theta
            Float theta = 0.5 * atan2(Sw(0, 1) + Sw(1, 0), Sw(0, 0) - Sw(1, 1));
            Float Ctheta = cos(theta);
            Float Stheta = sin(theta);
            Mat<2, 2> W;
            W << Ctheta, -Stheta, Stheta, Ctheta;

            // Compute singular values from Su
            Float SUsum = Su(0, 0) + Su(1, 1);
            Float SUdif = sqrt((Su(0, 0) - Su(1, 1)) * (Su(0, 0) - Su(1, 1)) + 4 * Su(0, 1) * Su(1, 0));
            Sigma(0) = sqrt((SUsum + SUdif) / 2);
            Sigma(1) = sqrt((SUsum - SUdif) / 2);

            // Find the correction matrix for the right side
            Mat<2, 2> S = U.transpose() * A * W;
            Eigen::DiagonalMatrix<Float, 2> C;
            C.diagonal() << ((S(0, 0) >= 0) ? 1 : -1), ((S(1, 1) >= 0) ? 1 : -1);
            V = W * C;
        }

        inline void svd_3x2(const Mat<3, 2>& A, Mat<3, 3>& U, Vec<2>& Sigma, Mat<2, 2>& V) {
            Mat<2, 2> U_;
            Mat<3, 2> proj;
            proj.col(0) = A.col(0).normalized();
            Vec3 normal = A.col(0).cross(A.col(1)).normalized();
            proj.col(1) = normal.cross(proj.col(0));
            Mat<2, 2> A_ = proj.transpose() * A;
            svd_2x2(A_, U_, Sigma, V);
            U.block(0, 0, 3, 2) = proj * U_;
            U.col(2) = normal;
        }

        inline Float invariant2(const Vec2& Sigma) {
            return Sigma[0] * Sigma[0] + Sigma[1] * Sigma[1];
        }

        inline void eigensystem_2x2(const Mat<2, 2>& A, Mat<2, 2>& Q, Vec2& Lambda) {
            // Find phi
            Float phi = 0.5 * atan2(A(0, 1) + A(1, 0), A(0, 0) - A(1, 1));
            Float Cphi = cos(phi);
            Float Sphi = sin(phi);
            Q << Cphi, -Sphi, Sphi, Cphi;
            Lambda(0) = A(0, 0) * Cphi * Cphi + A(1, 1) * Sphi * Sphi + (A(0, 1) + A(1, 0)) * Cphi * Sphi;
            Lambda(1) = A(0, 0) * Sphi * Sphi + A(1, 1) * Cphi * Cphi - (A(0, 1) + A(1, 0)) * Cphi * Sphi;
        }

        inline void buildScalingEigenvectors(const Mat<3, 3>& U, const Mat<2, 2>& Q, const Mat<2, 2>& V, Mat<6, 6>& Q6) {
            const Vec2 q0 = Q.col(0);
            const Vec2 q1 = Q.col(1);

            Mat<3, 2> q = Mat<3, 2>::Zero();
            q(0, 0) = q0(0);
            q(1, 1) = q0(1);
            const Mat<3, 2> Q0 = U * q * V.transpose();
            q(0, 0) = q1(0);
            q(1, 1) = q1(1);
            const Mat<3, 2> Q1 = U * q * V.transpose();

            Q6.col(4) = Eigen::Map<const Vec<6>>(Q0.data());
            Q6.col(5) = Eigen::Map<const Vec<6>>(Q1.data());
        }

        ///////////////////////////////////////////////////////////////////////
        // eigenvectors 0 is the twist mode
        // eigenvectors 1 is the flip mode
        ///////////////////////////////////////////////////////////////////////
        inline void buildTwistAndFlipEigenvectors(const Mat<3, 3>& U, const Mat<2, 2>& V, Mat<6, 6>& Q) {
            // create the twist matrices
            Mat<3, 2> T;
            T << 0, -1, 1, 0, 0, 0;

            const Mat<3, 2> Q0 = 0.70710678118 * (U * T * V.transpose());

            // create the flip matrices
            Mat<3, 2> L;
            L << 0, 1, 1, 0, 0, 0;

            const Mat<3, 2> Q1 = 0.70710678118 * (U * L * V.transpose());

            const Mat<3, 2> Q2 = U.col(2) * V.col(0).transpose();
            const Mat<3, 2> Q3 = U.col(2) * V.col(1).transpose();
            Q.col(0) = Eigen::Map<const Vec<6>>(Q0.data());
            Q.col(1) = Eigen::Map<const Vec<6>>(Q1.data());
            Q.col(2) = Eigen::Map<const Vec<6>>(Q2.data());
            Q.col(3) = Eigen::Map<const Vec<6>>(Q3.data());
        }

        inline void buildEigensystem(
            const Mat<3, 3>& U, const Vec2& Sigma, const Mat<2, 2>& V, Vec<6>& eigenvalues, Mat<6, 6>& eigenvectors, Float mu, Float lambda
        ) {
            const auto& _mu = mu;
            const auto& _lambda = lambda;
            const Float I2 = invariant2(Sigma);
            const Float front = _mu * (I2 - 1.0) + _lambda * 0.5 * (I2 - 2.0);
            const Float s0s1 = Sigma[0] * Sigma[1];
            const Float s0Sq = Sigma[0] * Sigma[0];
            const Float s1Sq = Sigma[1] * Sigma[1];

            // twists
            eigenvalues[0] = front - _mu * s0s1;

            // flips
            eigenvalues[1] = front + _mu * s0s1;

            eigenvalues[2] = front - _mu * s1Sq;
            eigenvalues[3] = front - _mu * s0Sq;

            // populate the scaling mode matrix
            Mat<2, 2> A;
            A(0, 0) = front + _lambda * s0Sq + _mu * (2 * s0Sq - s1Sq);
            A(1, 1) = front + _lambda * s1Sq + _mu * (2 * s1Sq - s0Sq);
            A(0, 1) = _lambda * s0s1;
            A(1, 0) = _lambda * s0s1;

            // get the scaling modes
            Mat<2, 2> Q;
            Vec2 Lambda;
            eigensystem_2x2(A, Q, Lambda);

            eigenvalues[4] = Lambda[0];
            eigenvalues[5] = Lambda[1];

            // Compute the eigenvectors
            buildTwistAndFlipEigenvectors(U, V, eigenvectors);
            buildScalingEigenvectors(U, Q, V, eigenvectors);
        }

        // 3. Second derivative (Hessian) w.r.t. F
        Mat<6, 6> PD2_StVKSsEnergy_DeformGrad(const Mat<3, 2> F, const SolverConfig& config)
        {
            Mat<6, 6> hessian = Mat<6, 6>::Zero();

            Float lambda = config.lambda;
            Float mu = config.mu;
            Float thickness = config.thickness;

            const Mat<2, 2> E = 0.5 * (F.transpose() * F - Mat<2, 2>::Identity());
            const auto traceE = E.trace();
            const Float psi = mu * E.squaredNorm() + 0.5 * lambda * traceE * traceE;

            const Mat<2, 2> dPsi_dE = 2.0 * mu * E + lambda * traceE * Mat<2, 2>::Identity();

            Mat<3, 3> U;
            Mat<2, 2> V;
            Vec2 Sigma;
            svd_3x2(F, U, Sigma, V);
            Vec<6> eigenvalues;
            Mat<6, 6> eigenvectors;

            buildEigensystem(U, Sigma, V, eigenvalues, eigenvectors, mu, lambda);
             //std::cout << eigenvalues << std::endl;
             //std::cout << eigenvectors << std::endl;
            for (int i = 0; i < 6; ++i) {
                if (eigenvalues[i] < 0)
                    eigenvalues[i] = 0;
            }
            hessian = eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
            hessian = hessian * thickness;

            //hessian = SPDProjection(hessian);

            return hessian;
        }

        Vec<27> StVKSsGradBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface)
        {
            Vec<27> mat = Vec<27>::Zero();       // support contains 9 nodes, each with 3 coordinate entries

            const Mat<3, 2>& F = quadPoint.deformationGradient;

            //this->AppendToLog(ToStr(F));

            DMat pd_SsEnergy_DeformGrad = PD_StVKSsEnergy_DeformGrad(F, config);

            DMat pd_DeformGrad_CtrlPt = PD_DeformGrad_CtrlPt_Block(quadPoint.U(), quadPoint.V(), indices, surface);

            mat = pd_DeformGrad_CtrlPt * pd_SsEnergy_DeformGrad;

            return config.thickness * mat;
        }
            
        Mat<27, 27> StVKSsHessBlockAt(const QuadPoint& quadPoint, const std::vector<UInt> indices, const SolverConfig& config, const BSSurface* const surface)
        {
            Mat<27, 27> mat = Mat<27, 27>::Zero();           // support contains 9 nodes, each with 3 coordinate entries

            const Mat<3, 2>& F = quadPoint.deformationGradient;

            DMat pd2_SsEnergy_DeformGrad = PD2_StVKSsEnergy_DeformGrad(F, config);
            DMat pd_DeformGrad_CtrlPt = PD_DeformGrad_CtrlPt_Block(quadPoint.U(), quadPoint.V(), indices, surface);

            // TODO test whether we can fill this in manually
            mat = pd_DeformGrad_CtrlPt * pd2_SsEnergy_DeformGrad * pd_DeformGrad_CtrlPt.transpose();
            Mat<27, 27> hess = config.thickness * mat;

            return hess;
        }
    }
}