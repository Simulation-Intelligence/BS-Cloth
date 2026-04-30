#pragma once

#pragma warning(push, 0)
#include <stdint.h>
#include <string>

#include "Eigen/Eigen"
#pragma warning(pop)

namespace BSIPC
{
    // STL Types
    using Int = int32_t;
    using UInt = uint32_t;
    using LInt = int64_t;
    using LUInt = uint64_t;
    using Bool = bool;
    using String = std::string;
    using Float = double;

    // Eigen Types
    using Vec2 = Eigen::Vector2d;
    using Vec3 = Eigen::Vector3d;
    using Vec4 = Eigen::Vector4d;

    using IVec2 = Eigen::Vector2i;
    using IVec3 = Eigen::Vector3i;
    using IVec4 = Eigen::Vector4i;

    using DVec = Eigen::VectorXd;
    
    template <Int X>
    using Vec = Eigen::Matrix<Float, X, 1, Eigen::ColMajor>;

    template <Int X, Int Y>
    using Mat = Eigen::Matrix<Float, X, Y, Eigen::ColMajor>;     // X = number of rows, Y = number of columns

    template <Int X, Int Y>
    using IMat = Eigen::Matrix<Int, X, Y, Eigen::ColMajor>;

    using DMat = Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

    using SpMat = Eigen::SparseMatrix<Float, Eigen::ColMajor, Int>;
    using SpMatData = std::vector<Eigen::Triplet<Float>>;
    using SpMatEntry = Eigen::Triplet<Float>;
}
