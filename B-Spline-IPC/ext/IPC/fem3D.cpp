#include "fem3D.h"
#include<iostream>
#include <vector>
#include <map>
#include "oneapi/tbb.h"
#include<fstream>
//#include <igl/winding_number.h>
#include "fem_parameters.h"
#include "IPCdistanceFuncs.h"

namespace IPC
{

using namespace std;
using namespace FEM;

namespace LibShell {

    Eigen::Matrix3d crossMatrix(Eigen::Vector3d v) {
        Eigen::Matrix3d ret;
        ret << 0, -v[2], v[1],
            v[2], 0, -v[0],
            -v[1], v[0], 0;
        return ret;
    }


    Eigen::Matrix2d adjugate(Eigen::Matrix2d M) {
        Eigen::Matrix2d ret;
        ret << M(1, 1), -M(0, 1), -M(1, 0), M(0, 0);
        return ret;
    }

    double angle(const Eigen::Vector3d& v, const Eigen::Vector3d& w, const Eigen::Vector3d& axis,
        Eigen::Matrix<double, 1, 9>* derivative, // v, w
        Eigen::Matrix<double, 9, 9>* hessian
    ) {
        double theta = 2.0 * atan2((v.cross(w).dot(axis) / axis.norm()), v.dot(w) + v.norm() * w.norm());

        if (derivative) {
            derivative->segment<3>(0) = -axis.cross(v) / v.squaredNorm() / axis.norm();
            derivative->segment<3>(3) = axis.cross(w) / w.squaredNorm() / axis.norm();
            derivative->segment<3>(6).setZero();
        }
        if (hessian) {
            hessian->setZero();
            hessian->block<3, 3>(0, 0) +=
                2.0 * (axis.cross(v)) * v.transpose() / v.squaredNorm() / v.squaredNorm() / axis.norm();
            hessian->block<3, 3>(3, 3) +=
                -2.0 * (axis.cross(w)) * w.transpose() / w.squaredNorm() / w.squaredNorm() / axis.norm();
            hessian->block<3, 3>(0, 0) += -crossMatrix(axis) / v.squaredNorm() / axis.norm();
            hessian->block<3, 3>(3, 3) += crossMatrix(axis) / w.squaredNorm() / axis.norm();

            Eigen::Matrix3d dahat = (Eigen::Matrix3d::Identity() / axis.norm() -
                axis * axis.transpose() / axis.norm() / axis.norm() / axis.norm());

            hessian->block<3, 3>(0, 6) += crossMatrix(v) * dahat / v.squaredNorm();
            hessian->block<3, 3>(3, 6) += -crossMatrix(w) * dahat / w.squaredNorm();
        }

        return theta;
    }
};

double edgeTheta(
    const Eigen::Vector3d& q0,
    const Eigen::Vector3d& q1,
    const Eigen::Vector3d& q2,
    const Eigen::Vector3d& q3,
    Eigen::Matrix<double, 1, 12>* derivative, // edgeVertex, then edgeOppositeVertex
    Eigen::Matrix<double, 12, 12>* hessian) {
    using namespace LibShell;
    if (derivative)
        derivative->setZero();
    if (hessian)
        hessian->setZero();

    Eigen::Vector3d n0 = (q0 - q2).cross(q1 - q2);
    Eigen::Vector3d n1 = (q1 - q3).cross(q0 - q3);
    Eigen::Vector3d axis = q1 - q0;
    Eigen::Matrix<double, 1, 9> angderiv;
    Eigen::Matrix<double, 9, 9> anghess;

    double theta = angle(n0, n1, axis, (derivative || hessian) ? &angderiv : NULL, hessian ? &anghess : NULL);

    if (derivative) {
        derivative->block<1, 3>(0, 0) += angderiv.block<1, 3>(0, 0) * crossMatrix(q2 - q1);
        derivative->block<1, 3>(0, 3) += angderiv.block<1, 3>(0, 0) * crossMatrix(q0 - q2);
        derivative->block<1, 3>(0, 6) += angderiv.block<1, 3>(0, 0) * crossMatrix(q1 - q0);

        derivative->block<1, 3>(0, 0) += angderiv.block<1, 3>(0, 3) * crossMatrix(q1 - q3);
        derivative->block<1, 3>(0, 3) += angderiv.block<1, 3>(0, 3) * crossMatrix(q3 - q0);
        derivative->block<1, 3>(0, 9) += angderiv.block<1, 3>(0, 3) * crossMatrix(q0 - q1);
    }

    if (hessian) {
        Eigen::Matrix3d vqm[3];
        vqm[0] = crossMatrix(q0 - q2);
        vqm[1] = crossMatrix(q1 - q0);
        vqm[2] = crossMatrix(q2 - q1);
        Eigen::Matrix3d wqm[3];
        wqm[0] = crossMatrix(q0 - q1);
        wqm[1] = crossMatrix(q1 - q3);
        wqm[2] = crossMatrix(q3 - q0);

        int vindices[3] = { 3, 6, 0 };
        int windices[3] = { 9, 0, 3 };

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                hessian->block<3, 3>(vindices[i], vindices[j]) +=
                    vqm[i].transpose() * anghess.block<3, 3>(0, 0) * vqm[j];
                hessian->block<3, 3>(vindices[i], windices[j]) +=
                    vqm[i].transpose() * anghess.block<3, 3>(0, 3) * wqm[j];
                hessian->block<3, 3>(windices[i], vindices[j]) +=
                    wqm[i].transpose() * anghess.block<3, 3>(3, 0) * vqm[j];
                hessian->block<3, 3>(windices[i], windices[j]) +=
                    wqm[i].transpose() * anghess.block<3, 3>(3, 3) * wqm[j];
            }

            hessian->block<3, 3>(vindices[i], 3) += vqm[i].transpose() * anghess.block<3, 3>(0, 6);
            hessian->block<3, 3>(3, vindices[i]) += anghess.block<3, 3>(6, 0) * vqm[i];
            hessian->block<3, 3>(vindices[i], 0) += -vqm[i].transpose() * anghess.block<3, 3>(0, 6);
            hessian->block<3, 3>(0, vindices[i]) += -anghess.block<3, 3>(6, 0) * vqm[i];

            hessian->block<3, 3>(windices[i], 3) += wqm[i].transpose() * anghess.block<3, 3>(3, 6);
            hessian->block<3, 3>(3, windices[i]) += anghess.block<3, 3>(6, 3) * wqm[i];
            hessian->block<3, 3>(windices[i], 0) += -wqm[i].transpose() * anghess.block<3, 3>(3, 6);
            hessian->block<3, 3>(0, windices[i]) += -anghess.block<3, 3>(6, 3) * wqm[i];

        }

        Eigen::Vector3d dang1 = angderiv.block<1, 3>(0, 0).transpose();
        Eigen::Vector3d dang2 = angderiv.block<1, 3>(0, 3).transpose();

        Eigen::Matrix3d dang1mat = crossMatrix(dang1);
        Eigen::Matrix3d dang2mat = crossMatrix(dang2);

        hessian->block<3, 3>(6, 3) += dang1mat;
        hessian->block<3, 3>(0, 3) -= dang1mat;
        hessian->block<3, 3>(0, 6) += dang1mat;
        hessian->block<3, 3>(3, 0) += dang1mat;
        hessian->block<3, 3>(3, 6) -= dang1mat;
        hessian->block<3, 3>(6, 0) -= dang1mat;

        hessian->block<3, 3>(9, 0) += dang2mat;
        hessian->block<3, 3>(3, 0) -= dang2mat;
        hessian->block<3, 3>(3, 9) += dang2mat;
        hessian->block<3, 3>(0, 3) += dang2mat;
        hessian->block<3, 3>(0, 9) -= dang2mat;
        hessian->block<3, 3>(9, 3) -= dang2mat;
    }

    return theta;
}






Vector4i calculateIsoDIndex_double(const vector<Vector3d> &vertexes, vector<uint64_t> &index) {
    int i = 0;
    int id1 = 1;
    int id2 = 2;
    int id3 = 3;

    Vector3d ed[3];

    ed[0] = vertexes[index[id1]] - vertexes[index[i]];
    ed[1] = vertexes[index[id2]] - vertexes[index[i]];
    ed[2] = vertexes[index[id3]] - vertexes[index[i]];

    double maxArea = 0;
    int mid = 0;
    for (int j = 0; j < 3; j++) {
        double area = ed[j].cross(ed[(j + 1) % 3]).norm();
        if (area > maxArea) {
            maxArea = area;
            mid = j;
        }
    }

    Vector3d ed0 = vertexes[index[id2]] - vertexes[index[id1]];
    Vector3d ed1 = vertexes[index[id3]] - vertexes[index[id1]];

    double area = ed0.cross(ed1).norm();
    if (area > maxArea) {
        maxArea = area;
        mid = 3;
    }
    int triangleId[3];
    if (mid < 3) {
        triangleId[0] = 0;
        triangleId[1] = mid + 1;
        triangleId[2] = (mid + 1) % 3 + 1;
        id2 = (mid + 2) % 3 + 1;
    } else {
        triangleId[0] = 1;
        triangleId[1] = 2;
        triangleId[2] = 3;
        id2 = 0;
    }
    double minL = std::numeric_limits<double>::max();
    int minId = 0;
    for (int j = 0; j < 3; j++) {
        Vector3d ted = vertexes[index[triangleId[(j + 1) % 3]]] - vertexes[index[triangleId[j]]];
        double L = ted.norm();
        if (L < minL) {
            minL = L;
            minId = j;
        }
    }
    i = triangleId[(minId + 2) % 3];
    id1 = triangleId[(minId + 1) % 3];
    id3 = triangleId[minId];


    //cout << i << " " << id1 << " " << id2 << " " << id3 << endl;
    ed[0] = vertexes[index[id1]] - vertexes[index[i]];
    ed[1] = vertexes[index[id2]] - vertexes[index[i]];
    ed[2] = vertexes[index[id3]] - vertexes[index[i]];

    Vector3d n = ed[2].cross(ed[0]);

    if (n.dot(ed[1]) < 0) {
        int temp = id1;
        id1 = id3;
        id3 = temp;
    }

    Vector4i temp_Index = Vector4i(index[i], index[id1], index[id2], index[id3]);
    index[0] = temp_Index[0];
    index[1] = temp_Index[1];
    index[2] = temp_Index[2];
    index[3] = temp_Index[3];
    Vector4i tet_id = Vector4i(i, id1, id2, id3);
    return tet_id;
}

void updateIndex_double(const vector<Vector3d> &vertexes, vector<uint64_t> &index) {
    int i = 0;
    int id1 = 1;
    int id2 = 2;
    int id3 = 3;

    Vector3d ed[3];

    ed[0] = vertexes[index[id1]] - vertexes[index[i]];
    ed[1] = vertexes[index[id2]] - vertexes[index[i]];
    ed[2] = vertexes[index[id3]] - vertexes[index[i]];

    double maxArea = 0;
    int mid = 0;
    for (int j = 0; j < 3; j++) {
        double area = ed[j].cross(ed[(j + 1) % 3]).norm();
        if (area > maxArea) {
            maxArea = area;
            mid = j;
        }
    }

    Vector3d ed0 = vertexes[index[id2]] - vertexes[index[id1]];
    Vector3d ed1 = vertexes[index[id3]] - vertexes[index[id1]];

    double area = ed0.cross(ed1).norm();
    if (area > maxArea) {
        maxArea = area;
        mid = 3;
    }
    int triangleId[3];
    if (mid < 3) {
        triangleId[0] = 0;
        triangleId[1] = mid + 1;
        triangleId[2] = (mid + 1) % 3 + 1;
        id2 = (mid + 2) % 3 + 1;
    } else {
        triangleId[0] = 1;
        triangleId[1] = 2;
        triangleId[2] = 3;
        id2 = 0;
    }
    double minL = std::numeric_limits<double>::max();
    int minId = 0;
    for (int j = 0; j < 3; j++) {
        Vector3d ted = vertexes[index[triangleId[(j + 1) % 3]]] - vertexes[index[triangleId[j]]];
        double L = ted.norm();
        if (L < minL) {
            minL = L;
            minId = j;
        }
    }
    i = triangleId[(minId + 2) % 3];
    id1 = triangleId[(minId + 1) % 3];
    id3 = triangleId[minId];


    //cout << i << " " << id1 << " " << id2 << " " << id3 << endl;
    ed[0] = vertexes[index[id1]] - vertexes[index[i]];
    ed[1] = vertexes[index[id2]] - vertexes[index[i]];
    ed[2] = vertexes[index[id3]] - vertexes[index[i]];

    Vector3d n = ed[2].cross(ed[0]);

    if (n.dot(ed[1]) < 0) {
        int temp = id1;
        id1 = id3;
        id3 = temp;
    }

    Vector4i temp_Index = Vector4i(index[i], index[id1], index[id2], index[id3]);
    index[0] = temp_Index[0];
    index[1] = temp_Index[1];
    index[2] = temp_Index[2];
    index[3] = temp_Index[3];
    //Vector4i tet_id = Vector4i(i, id1, id2, id3);
    //return tet_id;
}

void
solveLinearFunc(const Vector3d &e0, const Vector3d &e1, const Vector3d &e2, const Vector3d &N, double &b0, double &b2) {
    Matrix3d D, D1, D2;

    D << e0[0], e2[0], N[0], e0[1], e2[1], N[1], e0[2], e2[2], N[2];
    D1 << e1[0], e2[0], N[0], e1[1], e2[1], N[1], e1[2], e2[2], N[2];
    D2 << e0[0], e1[0], N[0], e0[1], e1[1], N[1], e0[2], e1[2], N[2];

    b0 = D1.determinant() / D.determinant();
    b2 = D2.determinant() / D.determinant();
}


Matrix3d
calculateIsoDms_double(const vector<Vector3d> &vertexes, const vector<uint64_t> &index, Vector3d &B, Vector3d *Normal,
                       bool isDs) {
    int i = 0;
    int id1 = 1;
    int id2 = 2;
    int id3 = 3;

    Vector3d ed[3];

    ed[0] = vertexes[index[id1]] - vertexes[index[i]];
    ed[1] = vertexes[index[id2]] - vertexes[index[i]];
    ed[2] = vertexes[index[id3]] - vertexes[index[i]];

    Vector3d n = (ed[2].cross(ed[0])).normalized();
    if (Normal != NULL) {
        *Normal = n;
    }
    double em;// = max(ed[0].norm(), ed[2].norm());

    double b0, b2;
    if (!isDs) {
        //Vector3d tb = ed[1] - ed[1].dot(n) * n;
        em = max(ed[0].norm(), ed[2].norm());
        solveLinearFunc(ed[0], ed[1], ed[2], n, b0, b2);
        /*if (abs(ed[0][0]) < 1e-15 && abs(ed[0][1]) < 1e-15) {
            if (abs(ed[0][2]) < 1e-15) {
                system("pause");
            }
            if (abs(tb[0]) > 1e-15) {
                b2 = tb[0] / ed[2][0];
                b0 = (tb[2] - b2 * ed[2][2]) / ed[0][2];
            }
            else if (abs(tb[1]) > 1e-15) {
                b2 = tb[1] / ed[2][1];
                b0 = (tb[2] - b2 * ed[2][2]) / ed[0][2];
            }
            else {
                b2 = 0;
                b0 = tb[2] / ed[0][2];
            }
        }
        else if (abs(ed[2][0]) < 1e-15 && abs(ed[2][1]) < 1e-15) {
            if (abs(ed[2][2]) < 1e-15) {
                system("pause");
            }
            if (abs(tb[0]) > 1e-15) {
                b0 = tb[0] / ed[0][0];
                b2 = (tb[2] - b0 * ed[0][2]) / ed[2][2];
            }
            else if (abs(tb[1]) > 1e-15) {
                b0 = tb[1] / ed[0][1];
                b2 = (tb[2] - b0 * ed[0][2]) / ed[2][2];
            }
            else {
                b0 = 0;
                b2 = tb[2] / ed[2][2];
            }
        }
        else if (abs(ed[0][0] * ed[2][1] - ed[2][0] * ed[0][1]) < 1e-15) {
            b2 = tb.dot(ed[2]) / (ed[2].norm()) / ed[2].norm();
            b0 = tb.dot(ed[0]) / (ed[0].norm()) / ed[0].norm();
        }
        else {

            b0 = (tb[0] * ed[2][1] - tb[1] * ed[2][0]) / (ed[0][0] * ed[2][1] - ed[2][0] * ed[0][1]);
            b2 = (tb[0] * ed[0][1] - tb[1] * ed[0][0]) / (ed[2][0] * ed[0][1] - ed[0][0] * ed[2][1]);
        }*/
        B[0] = b0;
        B[1] = b2;
        B[2] = em;
    } else {
        b0 = B[0];
        b2 = B[1];
        em = B[2];
    }


    Vector3d ned1 = b0 * ed[0] + b2 * ed[2] + n * em;

    Matrix3d M;
    M(0, 0) = ed[0][0];
    M(0, 1) = ned1[0];
    M(0, 2) = ed[2][0];
    M(1, 0) = ed[0][1];
    M(1, 1) = ned1[1];
    M(1, 2) = ed[2][1];
    M(2, 0) = ed[0][2];
    M(2, 1) = ned1[2];
    M(2, 2) = ed[2][2];

    return M;
}

void __Inverse2x2(const Eigen::Matrix2d &input, Eigen::Matrix2d &result) {
    double eps = 1e-15;
    const int dim = 2;
    double mat[dim * 2][dim * 2];
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < 2 * dim; j++) {
            if (j < dim) {
                mat[i][j] = input(i, j);
            } else {
                mat[i][j] = j - dim == i ? 1 : 0;
            }
        }
    }

    for (int i = 0; i < dim; i++) {
        if (abs(mat[i][i]) < eps) {
            int j;
            for (j = i + 1; j < dim; j++) {
                if (abs(mat[j][i]) > eps) break;
            }
            if (j == dim) return;
            for (int r = i; r < 2 * dim; r++) {
                mat[i][r] += mat[j][r];
            }
        }
        double ep = mat[i][i];
        for (int r = i; r < 2 * dim; r++) {
            mat[i][r] /= ep;
        }

        for (int j = i + 1; j < dim; j++) {
            double e = -1 * (mat[j][i] / mat[i][i]);
            for (int r = i; r < 2 * dim; r++) {
                mat[j][r] += e * mat[i][r];
            }
        }
    }

    for (int i = dim - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            double e = -1 * (mat[j][i] / mat[i][i]);
            for (int r = i; r < 2 * dim; r++) {
                mat[j][r] += e * mat[i][r];
            }
        }
    }


    for (int i = 0; i < dim; i++) {
        for (int r = dim; r < 2 * dim; r++) {
            result(i, r - dim) = mat[i][r];
        }
    }
}



void __Inverse(const Eigen::Matrix3d &input, Eigen::Matrix3d &result) {
    double eps = 1e-15;
    const int dim = 3;
    double mat[dim * 2][dim * 2];
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < 2 * dim; j++) {
            if (j < dim) {
                mat[i][j] = input(i, j);
            } else {
                mat[i][j] = j - dim == i ? 1 : 0;
            }
        }
    }

    for (int i = 0; i < dim; i++) {
        if (abs(mat[i][i]) < eps) {
            int j;
            for (j = i + 1; j < dim; j++) {
                if (abs(mat[j][i]) > eps) break;
            }
            if (j == dim) return;
            for (int r = i; r < 2 * dim; r++) {
                mat[i][r] += mat[j][r];
            }
        }
        double ep = mat[i][i];
        for (int r = i; r < 2 * dim; r++) {
            mat[i][r] /= ep;
        }

        for (int j = i + 1; j < dim; j++) {
            double e = -1 * (mat[j][i] / mat[i][i]);
            for (int r = i; r < 2 * dim; r++) {
                mat[j][r] += e * mat[i][r];
            }
        }
    }

    for (int i = dim - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            double e = -1 * (mat[j][i] / mat[i][i]);
            for (int r = i; r < 2 * dim; r++) {
                mat[j][r] += e * mat[i][r];
            }
        }
    }


    for (int i = 0; i < dim; i++) {
        for (int r = dim; r < 2 * dim; r++) {
            result(i, r - dim) = mat[i][r];
        }
    }
}


double calculateVolum(const vector<Vector3d> &vertexes, const Vector4i &index) {
    int id0 = 0;
    int id1 = 1;
    int id2 = 2;
    int id3 = 3;
    double o1x = vertexes[index[id1]][0] - vertexes[index[id0]][0];
    double o1y = vertexes[index[id1]][1] - vertexes[index[id0]][1];
    double o1z = vertexes[index[id1]][2] - vertexes[index[id0]][2];
    Vector3d OA = Vector3d(o1x, o1y, o1z);

    double o2x = vertexes[index[id2]][0] - vertexes[index[id0]][0];
    double o2y = vertexes[index[id2]][1] - vertexes[index[id0]][1];
    double o2z = vertexes[index[id2]][2] - vertexes[index[id0]][2];
    Vector3d OB = Vector3d(o2x, o2y, o2z);

    double o3x = vertexes[index[id3]][0] - vertexes[index[id0]][0];
    double o3y = vertexes[index[id3]][1] - vertexes[index[id0]][1];
    double o3z = vertexes[index[id3]][2] - vertexes[index[id0]][2];
    Vector3d OC = Vector3d(o3x, o3y, o3z);

    Vector3d heightDir = OA.cross(OB);
    double bottomArea = heightDir.norm();
    heightDir.normalize();

    double volum = bottomArea * abs(heightDir.dot(OC)) / 6;
    return volum > 0 ? volum : -volum;
}

double calculateTriangleArea(const vector<Vector3d> &vertexes, const Vector3i &index) {
    auto v01 = vertexes[index[1]] - vertexes[index[0]];
    auto v02 = vertexes[index[2]] - vertexes[index[0]];
    auto area = v01.cross(v02).norm() *0.5;
    return area;
}

Matrix3d calculateDms3D_double(const vector<Vector3d> &vertexes, const Vector4i &index, const int &i) {
    int id1 = (i + 1) % 4;
    int id2 = (i + 2) % 4;
    int id3 = (i + 3) % 4;
    double o1x = vertexes[index[id1]][0] - vertexes[index[i]][0];
    double o1y = vertexes[index[id1]][1] - vertexes[index[i]][1];
    double o1z = vertexes[index[id1]][2] - vertexes[index[i]][2];

    double o2x = vertexes[index[id2]][0] - vertexes[index[i]][0];
    double o2y = vertexes[index[id2]][1] - vertexes[index[i]][1];
    double o2z = vertexes[index[id2]][2] - vertexes[index[i]][2];

    double o3x = vertexes[index[id3]][0] - vertexes[index[i]][0];
    double o3y = vertexes[index[id3]][1] - vertexes[index[i]][1];
    double o3z = vertexes[index[id3]][2] - vertexes[index[i]][2];

    Matrix3d M;
    M(0, 0) = o1x;
    M(0, 1) = o2x;
    M(0, 2) = o3x;
    M(1, 0) = o1y;
    M(1, 1) = o2y;
    M(1, 2) = o3y;
    M(2, 0) = o1z;
    M(2, 1) = o2z;
    M(2, 2) = o3z;

    return M;
}


Matrix<double, 3, 2> calculateDs32D_double(const vector<Vector3d> &vertexes, const Vector3i &index) {
    Matrix<double, 3, 2> M;
    M.col(0) = vertexes[index[1]] - vertexes[index[0]];
    M.col(1) = vertexes[index[2]] - vertexes[index[0]];

    return M;
}


Matrix2d calculateDms2D_double(const vector<Vector3d> &vertexes, const Vector3i &index) {
    Eigen::Vector3d v01 = vertexes[index[1]] - vertexes[index[0]];
    Eigen::Vector3d v02 = vertexes[index[2]] - vertexes[index[0]];
    // compute uv coordinates by rotating each triangle normal to (0, 1, 0)
    Eigen::Vector3d normal = v01.cross(v02).normalized();
    Eigen::Vector3d target = Eigen::Vector3d(0, 1, 0);


    Eigen::Vector3d vec = normal.cross(target);
    double sin = vec.norm();
    double cos = normal.dot(target);
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d cross_vec;
    cross_vec << 0, -vec.z(), vec.y(),
            vec.z(), 0, -vec.x(),
            -vec.y(), vec.x(), 0;

    rotation += cross_vec + cross_vec * cross_vec / (1 + cos);

    auto rotate_uv0 = rotation * vertexes[index[0]];
    auto rotate_uv1 = rotation * vertexes[index[1]];
    auto rotate_uv2 = rotation * vertexes[index[2]];

    auto uv0 = Eigen::Vector2d(rotate_uv0.x(), rotate_uv0.z());
    auto uv1 = Eigen::Vector2d(rotate_uv1.x(), rotate_uv1.z());
    auto uv2 = Eigen::Vector2d(rotate_uv2.x(), rotate_uv2.z());
    Matrix2d M;
    M.col(0) = uv1 - uv0;
    M.col(1) = uv2 - uv0;
    return M;
}

Matrix3d calculateDms3D2_double(const vector<Vector3d> &vertexes, const vector<uint64_t> &index, const Vector4i &id) {
    int i = id[0];
    int id1 = id[1];
    int id2 = id[2];
    int id3 = id[3];
    double o1x = vertexes[index[id1]][0] - vertexes[index[i]][0];
    double o1y = vertexes[index[id1]][1] - vertexes[index[i]][1];
    double o1z = vertexes[index[id1]][2] - vertexes[index[i]][2];

    double o2x = vertexes[index[id2]][0] - vertexes[index[i]][0];
    double o2y = vertexes[index[id2]][1] - vertexes[index[i]][1];
    double o2z = vertexes[index[id2]][2] - vertexes[index[i]][2];

    double o3x = vertexes[index[id3]][0] - vertexes[index[i]][0];
    double o3y = vertexes[index[id3]][1] - vertexes[index[i]][1];
    double o3z = vertexes[index[id3]][2] - vertexes[index[i]][2];

    Matrix3d M;
    M(0, 0) = o1x;
    M(0, 1) = o2x;
    M(0, 2) = o3x;
    M(1, 0) = o1y;
    M(1, 1) = o2y;
    M(1, 2) = o3y;
    M(2, 0) = o1z;
    M(2, 1) = o2z;
    M(2, 2) = o3z;

    return M;
}

double calculateRestMaxEdg_double(const vector<Vector3d> &vertexes, const vector<uint64_t> &index, const int &i) {
    int id1 = (i + 1) % 4;
    int id2 = (i + 2) % 4;
    int id3 = (i + 3) % 4;
    double o1x = vertexes[index[id1]][0] - vertexes[index[i]][0];
    double o1y = vertexes[index[id1]][1] - vertexes[index[i]][1];
    double o1z = vertexes[index[id1]][2] - vertexes[index[i]][2];
    Vector3d ed0 = Vector3d(o1x, o1y, o1z);

    double o2x = vertexes[index[id2]][0] - vertexes[index[i]][0];
    double o2y = vertexes[index[id2]][1] - vertexes[index[i]][1];
    double o2z = vertexes[index[id2]][2] - vertexes[index[i]][2];
    Vector3d ed1 = Vector3d(o2x, o2y, o2z);

    double o3x = vertexes[index[id3]][0] - vertexes[index[i]][0];
    double o3y = vertexes[index[id3]][1] - vertexes[index[i]][1];
    double o3z = vertexes[index[id3]][2] - vertexes[index[i]][2];
    Vector3d ed2 = Vector3d(o3x, o3y, o3z);

    return max(max(ed0.norm(), ed1.norm()), ed2.norm());
}

//double getObjEnergy_StableNHK2_3D(const vector<Vector3d> &vertexes, const mesh3D &mesh, const double &lengthRate,
//                                  const double &volumeRate) {
//    double energy = parallel_reduce(
//            tbb::blocked_range<int>(0, mesh.tetrahedraNum), 0.0,
//            [&](const tbb::blocked_range<int> &rg, double temp_sum) {
//                for (int i = rg.begin(); i != rg.end(); i++) {
//                    MatrixXd F =
//                            calculateDms3D_double(vertexes, mesh.tetrahedras[i], 0) * mesh.DM_tetrahedra_inverse[i];
//
//                    SVDResult3D_double svdResult = QRSVD(F);
//                    Matrix3d V, S, sigma;
//                    //U = svdResult.U;
//                    sigma = svdResult.SIGMA;
//
//
//                    V = svdResult.V;
//
//                    S = V * sigma * V.transpose();
//                    double I2 = (S * S).trace();
//                    double I3 = S.determinant();
//                    temp_sum += (0.5 * lengthRate * (I2 - 3) +
//                                 0.5 * volumeRate * (I3 - 1 - 3 * lengthRate / 4 / volumeRate) *
//                                 (I3 - 1 - 3 * lengthRate / 4 / volumeRate) - 0.5 * lengthRate * log(I2 + 1)) *
//                                mesh.volum[i];
//
//                }
//                return temp_sum;
//            },
//            [&](double left, double right) {
//                return left + right;
//            }
//    );
//    return energy;
//}

double getObjBending_Energy(const mesh3D& mesh) {
    double energy = parallel_reduce(
        tbb::blocked_range<int>(0, mesh.tri_edges.size()), 0.0,
        [&](const tbb::blocked_range<int>& rg, double temp_sum) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                
                double energyVal = 0;
                const auto& adj = mesh.tri_edges_adj_points[i];
                const auto& edge = mesh.tri_edges[i];
                if (adj.y() != -1) {

                    auto x0 = mesh.vertices[edge.x()];
                    auto x1 = mesh.vertices[edge.y()];
                    auto x2 = mesh.vertices[adj.x()];
                    auto x3 = mesh.vertices[adj.y()];

                    double t = edgeTheta(x0, x1, x2, x3, nullptr, nullptr);

                    auto rest_x0 = mesh.v_rest[edge.x()];
                    auto rest_x1 = mesh.v_rest[edge.y()];
                    auto rest_x2 = mesh.v_rest[adj.x()];
                    auto rest_x3 = mesh.v_rest[adj.y()];

                    double length = (rest_x0 - rest_x1).norm();

                    double rest_t = edgeTheta(rest_x0, rest_x1, rest_x2, rest_x3, nullptr, nullptr);
                    energyVal = mesh.bendingStiffness * (t - rest_t) * (t - rest_t) * length;
                    
                }
                temp_sum += energyVal;

            }
            return temp_sum;
        },
        [&](double left, double right) {
            return left + right;
        }
        );
    return energy;
}

double getObjEnergy_baraffwitkin_3D(const mesh3D& mesh,
    const Vector2d& anisotropic_a,
    const Vector2d& anisotropic_b,
    double stretchS, double shearS, double strainRate) {
    double energy = parallel_reduce(
        tbb::blocked_range<int>(0, mesh.triangleNum), 0.0,
        [&](const tbb::blocked_range<int>& rg, double temp_sum) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                Matrix<double, 3, 2> F =
                    calculateDs32D_double(mesh.vertices, mesh.triangles[i]) * mesh.DM_triangle_inverse[i];
                double I6 = anisotropic_a.transpose() * F.transpose() * F * anisotropic_b;
                double shear_energy = I6 * I6;


                double I5u = (F * anisotropic_a).norm();
                double I5v = (F * anisotropic_b).norm();

                double ucoeff = 1;
                double vcoeff = 1;

                if (I5u <= 1) {
                    ucoeff = 0;
                }
                if (I5v <= 1) {
                    vcoeff = 0;
                }



                double stretch_energy =
                    pow(I5u - 1, 2) + ucoeff * strainRate * pow(I5u - 1, 3) +
                    pow(I5v - 1, 2) + vcoeff * strainRate * pow(I5v - 1, 3);
                temp_sum += (stretchS * stretch_energy + shearS * shear_energy) * mesh.areas[i];

            }
            return temp_sum;
        },
        [&](double left, double right) {
            return left + right;
        }
    );
    return energy;
}

//double getObjRestEnergy_StableNHK2_3D(const vector<Vector3d> &vertexes, const mesh3D &mesh, const double &lengthRate,
//                                      const double &volumeRate) {
//    double energy = parallel_reduce(
//            tbb::blocked_range<int>(0, mesh.tetrahedraNum), 0.0,
//            [&](const tbb::blocked_range<int> &rg, double temp_sum) {
//                for (int i = rg.begin(); i != rg.end(); i++) {
//                    MatrixXd F =
//                            calculateDms3D_double(vertexes, mesh.tetrahedras[i], 0) * mesh.DM_tetrahedra_inverse[i];
//                    SVDResult3D_double svdResult = QRSVD(F);
//                    Matrix3d V, S, sigma;
//                    sigma = svdResult.SIGMA;
//                    V = svdResult.V;
//
//                    S = V * sigma * V.transpose();
//                    double I2 = (S * S).trace();
//                    double I3 = S.determinant();
//                    temp_sum +=
//                            (0.5 * volumeRate * (3 * lengthRate / 4 / volumeRate) * (3 * lengthRate / 4 / volumeRate) -
//                             0.5 * lengthRate * log(4)) * mesh.volum[i];
//                }
//                return temp_sum;
//            },
//            [&](double left, double right) {
//                return left + right;
//            }
//    );
//    return energy;
//}


MatrixXd computePFPX3D_double(const Matrix3d &InverseDm) {
    MatrixXd PFPX = MatrixXd::Zero(9, 12);
    double m = InverseDm(0, 0), n = InverseDm(0, 1), o = InverseDm(0, 2);
    double p = InverseDm(1, 0), q = InverseDm(1, 1), r = InverseDm(1, 2);
    double s = InverseDm(2, 0), t = InverseDm(2, 1), u = InverseDm(2, 2);
    double t1 = -(m + p + s);
    double t2 = -(n + q + t);
    double t3 = -(o + r + u);
    PFPX(0, 0) = t1;
    PFPX(0, 3) = m;
    PFPX(0, 6) = p;
    PFPX(0, 9) = s;
    PFPX(1, 1) = t1;
    PFPX(1, 4) = m;
    PFPX(1, 7) = p;
    PFPX(1, 10) = s;
    PFPX(2, 2) = t1;
    PFPX(2, 5) = m;
    PFPX(2, 8) = p;
    PFPX(2, 11) = s;
    PFPX(3, 0) = t2;
    PFPX(3, 3) = n;
    PFPX(3, 6) = q;
    PFPX(3, 9) = t;
    PFPX(4, 1) = t2;
    PFPX(4, 4) = n;
    PFPX(4, 7) = q;
    PFPX(4, 10) = t;
    PFPX(5, 2) = t2;
    PFPX(5, 5) = n;
    PFPX(5, 8) = q;
    PFPX(5, 11) = t;
    PFPX(6, 0) = t3;
    PFPX(6, 3) = o;
    PFPX(6, 6) = r;
    PFPX(6, 9) = u;
    PFPX(7, 1) = t3;
    PFPX(7, 4) = o;
    PFPX(7, 7) = r;
    PFPX(7, 10) = u;
    PFPX(8, 2) = t3;
    PFPX(8, 5) = o;
    PFPX(8, 8) = r;
    PFPX(8, 11) = u;
    return PFPX;
}


MatrixXd computePFPX32D_double(const Matrix2d &InverseDm) {
    double s0 = InverseDm.col(0).sum();
    double s1 = InverseDm.col(1).sum();

    double d0 = InverseDm(0, 0);
    double d1 = InverseDm(1, 0);
    double d2 = InverseDm(0, 1);
    double d3 = InverseDm(1, 1);
    Eigen::Matrix<double, 6, 9> dFdx;
    dFdx.setZero();
    // dF / dx0
    dFdx(0, 0) = -s0;
    dFdx(3, 0) = -s1;

    // dF / dy0
    dFdx(1, 1) = -s0;
    dFdx(4, 1) = -s1;

    // dF / dz0
    dFdx(2, 2) = -s0;
    dFdx(5, 2) = -s1;

    // dF / dx1
    dFdx(0, 3) = d0;
    dFdx(3, 3) = d2;

    // dF / dy1
    dFdx(1, 4) = d0;
    dFdx(4, 4) = d2;

    // dF / dz1
    dFdx(2, 5) = d0;
    dFdx(5, 5) = d2;

    // dF / dx2
    dFdx(0, 6) = d1;
    dFdx(3, 6) = d3;

    // dF / dy2
    dFdx(1, 7) = d1;
    dFdx(4, 7) = d3;

    // dF / dz2
    dFdx(2, 8) = d1;
    dFdx(5, 8) = d3;
    return dFdx;
}


MatrixXd computeIsoPFPX3D_double(const Matrix3d &InverseDm, const Vector3d &B, const vector<Vector3d> &vertexes,
                                 const vector<uint64_t> &index) {
    MatrixXd PFPX = MatrixXd::Zero(9, 12);
    int i = 0;
    int id1 = 1;
    int id2 = 2;
    int id3 = 3;
    double b0 = B[0], b2 = B[1], emax = B[2];
    Vector3d ed[3];

    ed[0] = vertexes[index[id1]] - vertexes[index[i]];
    ed[1] = vertexes[index[id2]] - vertexes[index[i]];
    ed[2] = vertexes[index[id3]] - vertexes[index[i]];

    //Vector3d N = (ed[2].cross(ed[0])).normalized();

    Vector3d dndu[12];
    dndu[0] = Vector3d(0, ed[0][2] - ed[2][2], ed[2][1] - ed[0][1]);
    dndu[1] = Vector3d(ed[2][2] - ed[0][2], 0, ed[0][0] - ed[2][0]);
    dndu[2] = Vector3d(ed[0][1] - ed[2][1], ed[2][0] - ed[0][0], 0);
    dndu[3] = Vector3d(0, ed[2][2], -ed[2][1]);
    dndu[4] = Vector3d(-ed[2][2], 0, ed[2][0]);
    dndu[5] = Vector3d(ed[2][1], -ed[2][0], 0);
    dndu[6] = Vector3d(0, 0, 0);
    dndu[7] = Vector3d(0, 0, 0);
    dndu[8] = Vector3d(0, 0, 0);
    dndu[9] = Vector3d(0, -ed[0][2], ed[0][1]);
    dndu[10] = Vector3d(ed[0][2], 0, -ed[0][0]);
    dndu[11] = Vector3d(-ed[0][1], ed[0][0], 0);

    for (int i = 0; i < 12; i++) {
        Vector3d ed20 = ed[2].cross(ed[0]);
        double tempR = ed20.transpose() * dndu[i];
        //cout << ed20.transpose() * dndu[i] << endl;
        double tempR2 = ed20.transpose() * ed20;
        dndu[i] = (1.0 / (ed20.norm()) * dndu[i] - tempR / pow(tempR2, 1.5) * ed20);
        dndu[i] = emax * dndu[i];
    }

    double m = InverseDm(0, 0), n = InverseDm(0, 1), o = InverseDm(0, 2);
    double p = InverseDm(1, 0), q = InverseDm(1, 1), r = InverseDm(1, 2);
    double s = InverseDm(2, 0), t = InverseDm(2, 1), u = InverseDm(2, 2);


    PFPX(0, 0) = -m - s + (dndu[0][0] - b0 - b2) * p;
    PFPX(3, 0) = -n - t + (dndu[0][0] - b0 - b2) * q;
    PFPX(6, 0) = -o - u + (dndu[0][0] - b0 - b2) * r;
    PFPX(1, 0) = p * dndu[0][1];
    PFPX(4, 0) = q * dndu[0][1];
    PFPX(7, 0) = r * dndu[0][1];
    PFPX(2, 0) = p * dndu[0][2];
    PFPX(5, 0) = q * dndu[0][2];
    PFPX(8, 0) = r * dndu[0][2];

    PFPX(0, 1) = dndu[1][0] * p;
    PFPX(3, 1) = dndu[1][0] * q;
    PFPX(6, 1) = dndu[1][0] * r;
    PFPX(1, 1) = -m - s + (dndu[1][1] - b0 - b2) * p;
    PFPX(4, 1) = -n - t + (dndu[1][1] - b0 - b2) * q;
    PFPX(7, 1) = -o - u + (dndu[1][1] - b0 - b2) * r;
    PFPX(2, 1) = dndu[1][2] * p;
    PFPX(5, 1) = dndu[1][2] * q;
    PFPX(8, 1) = dndu[1][2] * r;

    PFPX(0, 2) = dndu[2][0] * p;
    PFPX(3, 2) = dndu[2][0] * q;
    PFPX(6, 2) = dndu[2][0] * r;
    PFPX(1, 2) = dndu[2][1] * p;
    PFPX(4, 2) = dndu[2][1] * q;
    PFPX(7, 2) = dndu[2][1] * r;
    PFPX(2, 2) = -m - s + (dndu[2][2] - b0 - b2) * p;
    PFPX(5, 2) = -n - t + (dndu[2][2] - b0 - b2) * q;
    PFPX(8, 2) = -o - u + (dndu[2][2] - b0 - b2) * r;


    PFPX(0, 3) = m + (dndu[3][0] + b0) * p;
    PFPX(3, 3) = n + (dndu[3][0] + b0) * q;
    PFPX(6, 3) = o + (dndu[3][0] + b0) * r;
    PFPX(1, 3) = p * dndu[3][1];
    PFPX(4, 3) = q * dndu[3][1];
    PFPX(7, 3) = r * dndu[3][1];
    PFPX(2, 3) = p * dndu[3][2];
    PFPX(5, 3) = q * dndu[3][2];
    PFPX(8, 3) = r * dndu[3][2];

    PFPX(0, 4) = dndu[4][0] * p;
    PFPX(3, 4) = dndu[4][0] * q;
    PFPX(6, 4) = dndu[4][0] * r;
    PFPX(1, 4) = m + (dndu[4][1] + b0) * p;
    PFPX(4, 4) = n + (dndu[4][1] + b0) * q;
    PFPX(7, 4) = o + (dndu[4][1] + b0) * r;
    PFPX(2, 4) = dndu[4][2] * p;
    PFPX(5, 4) = dndu[4][2] * q;
    PFPX(8, 4) = dndu[4][2] * r;

    PFPX(0, 5) = dndu[5][0] * p;
    PFPX(3, 5) = dndu[5][0] * q;
    PFPX(6, 5) = dndu[5][0] * r;
    PFPX(1, 5) = dndu[5][1] * p;
    PFPX(4, 5) = dndu[5][1] * q;
    PFPX(7, 5) = dndu[5][1] * r;
    PFPX(2, 5) = m + (dndu[5][2] + b0) * p;
    PFPX(5, 5) = n + (dndu[5][2] + b0) * q;
    PFPX(8, 5) = o + (dndu[5][2] + b0) * r;


    PFPX(0, 9) = s + (dndu[9][0] + b2) * p;
    PFPX(3, 9) = t + (dndu[9][0] + b2) * q;
    PFPX(6, 9) = u + (dndu[9][0] + b2) * r;
    PFPX(1, 9) = p * dndu[9][1];
    PFPX(4, 9) = q * dndu[9][1];
    PFPX(7, 9) = r * dndu[9][1];
    PFPX(2, 9) = p * dndu[9][2];
    PFPX(5, 9) = q * dndu[9][2];
    PFPX(8, 9) = r * dndu[9][2];

    PFPX(0, 10) = dndu[10][0] * p;
    PFPX(3, 10) = dndu[10][0] * q;
    PFPX(6, 10) = dndu[10][0] * r;
    PFPX(1, 10) = s + (dndu[10][1] + b2) * p;
    PFPX(4, 10) = t + (dndu[10][1] + b2) * q;
    PFPX(7, 10) = u + (dndu[10][1] + b2) * r;
    PFPX(2, 10) = dndu[10][2] * p;
    PFPX(5, 10) = dndu[10][2] * q;
    PFPX(8, 10) = dndu[10][2] * r;

    PFPX(0, 11) = dndu[11][0] * p;
    PFPX(3, 11) = dndu[11][0] * q;
    PFPX(6, 11) = dndu[11][0] * r;
    PFPX(1, 11) = dndu[11][1] * p;
    PFPX(4, 11) = dndu[11][1] * q;
    PFPX(7, 11) = dndu[11][1] * r;
    PFPX(2, 11) = s + (dndu[11][2] + b2) * p;
    PFPX(5, 11) = t + (dndu[11][2] + b2) * q;
    PFPX(8, 11) = u + (dndu[11][2] + b2) * r;

    return PFPX;
}


double getRandValue(const double &left, const double &right) {
    double value = (right - left) * rand() / double(RAND_MAX) + left;
    return value;
}

float distance(Vector3d x, Vector3d y) {
    return (x - y).norm();
}


void initMesh3D(mesh3D &mesh, int type, double scale) {

    mesh.masses = vector<double>(mesh.vertexNum, 0);
    //double min = 10000000000000;
    double massSum = 0;
    for (int i = 0; i < mesh.tetrahedralNum; i++) {
        Matrix3d DM = calculateDms3D_double(mesh.vertices, mesh.tetrahedrals[i], 0);
        Matrix3d ONE;
        //ONE << 1, 1, 1, 1, 1, 1, 1, 1, 1;
        double vlm = calculateVolum(mesh.vertices, mesh.tetrahedrals[i]);


        mesh.masses[mesh.tetrahedrals[i][0]] += vlm * mesh.density / 4;
        mesh.masses[mesh.tetrahedrals[i][1]] += vlm * mesh.density / 4;
        mesh.masses[mesh.tetrahedrals[i][2]] += vlm * mesh.density / 4;
        mesh.masses[mesh.tetrahedrals[i][3]] += vlm * mesh.density / 4;

        massSum += vlm * mesh.density;

        Eigen::Matrix3d DMInverse = DM.inverse();
//        __Inverse(DM, DMInverse);
//        mesh.DM_tetrahedra_inverse.push_back(DMInverse);

        mesh.DM_tetrahedra_inverse.push_back(DMInverse);
        mesh.volum.push_back(vlm);

    }
    // triangles
    for (int i = 0; i < mesh.triangleNum; i++) {
        Matrix2d DM = calculateDms2D_double(mesh.vertices, mesh.triangles[i]);
//        Matrix3d ONE;
        //ONE << 1, 1, 1, 1, 1, 1, 1, 1, 1;
        double area = calculateTriangleArea(mesh.vertices, mesh.triangles[i]);
        
        area *= mesh.clothThickness;

        mesh.areas.push_back(area);

        mesh.masses[mesh.triangles[i][0]] += area * mesh.cloth_density / 3;
        mesh.masses[mesh.triangles[i][1]] += area * mesh.cloth_density / 3;
        mesh.masses[mesh.triangles[i][2]] += area * mesh.cloth_density / 3;

        massSum += area * mesh.cloth_density;


//        Eigen::Matrix2d DMInverse;
//        __Inverse2x2(DM, DMInverse);
//        mesh.DM_triangle_inverse.push_back(DMInverse);
//
        Eigen::Matrix2d DMInverse = DM.inverse();
        mesh.DM_triangle_inverse.push_back(DMInverse);

        //mesh.DM_tetrahedra_inverse.push_back(DM.inverse());


    }

    mesh.meanMass = massSum / mesh.vertexNum;
    printf("meanMass: %f\n", mesh.meanMass);
}


//Matrix3d computePEPF_StableNHK3D_double(const Matrix3d &F, const double &lengthRate, const double &volumRate) {
//    SVDResult3D_double svdResult = QRSVD(F);
//    Matrix3d U, V, R, S, sigma;
//    //U = svdResult.U;
//    sigma = svdResult.SIGMA;
//    V = svdResult.V;
//
//    S = V * sigma * V.transpose();
//
//    double u = lengthRate, r = volumRate;
//    Matrix3d pI3pF;
//    pI3pF(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
//    pI3pF(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
//    pI3pF(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);
//
//    pI3pF(1, 0) = F(2, 1) * F(0, 2) - F(2, 2) * F(0, 1);
//    pI3pF(1, 1) = F(2, 2) * F(0, 0) - F(2, 0) * F(0, 2);
//    pI3pF(1, 2) = F(2, 0) * F(0, 1) - F(2, 1) * F(0, 0);
//
//    pI3pF(2, 0) = F(0, 1) * F(1, 2) - F(1, 1) * F(0, 2);
//    pI3pF(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
//    pI3pF(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
//    Matrix3d PEPF = u * F + (r * (S.determinant() - 1) - u) * pI3pF;
//    return PEPF;
//}
//
//Matrix3d computePEPF_StableNHK3D_2_double(const Matrix3d &F, const double &lengthRate, const double &volumRate) {
//    SVDResult3D_double svdResult = QRSVD(F);
//    Matrix3d U, V, R, S, sigma;
//
//    sigma = svdResult.SIGMA;
//
//    double I3 = sigma(0, 0) * sigma(1, 1) * sigma(2, 2);
//    double I2 = sigma(0, 0) * sigma(0, 0) + sigma(1, 1) * sigma(1, 1) + sigma(2, 2) * sigma(2, 2);
//
//    double u = lengthRate, r = volumRate;
//    Matrix3d pI3pF;
//
//    pI3pF(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
//    pI3pF(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
//    pI3pF(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);
//
//    pI3pF(1, 0) = F(2, 1) * F(0, 2) - F(2, 2) * F(0, 1);
//    pI3pF(1, 1) = F(2, 2) * F(0, 0) - F(2, 0) * F(0, 2);
//    pI3pF(1, 2) = F(2, 0) * F(0, 1) - F(2, 1) * F(0, 0);
//
//    pI3pF(2, 0) = F(0, 1) * F(1, 2) - F(1, 1) * F(0, 2);
//    pI3pF(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
//    pI3pF(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
//
//    Matrix3d PEPF = u * (1 - 1 / (I2 + 1)) * F + (r * (I3 - 1 - u * 3 / (r * 4))) * pI3pF;
//    return PEPF;
//}


Matrix<double, 3, 2>
computePEPF_baraffwitkin_double(const Matrix<double, 3, 2>& F,
    const Vector2d& anisotropic_a,
    const Vector2d& anisotropic_b,
    double stretchS, double shearS, double strainRate) {
    double I6 = anisotropic_a.transpose() * F.transpose() * F * anisotropic_b;
    Eigen::Matrix<double, 3, 2> stretch_pk1, shear_pk1;

    shear_pk1 = 2 * (I6 - anisotropic_a.transpose() * anisotropic_b) *
        (F * anisotropic_a * anisotropic_b.transpose() +
            F * anisotropic_b * anisotropic_a.transpose()
            );
    double I5u = (F * anisotropic_a).transpose() * F * anisotropic_a;
    double I5v = (F * anisotropic_b).transpose() * F * anisotropic_b;
    double ucoeff = 1.0 - 1 / sqrt(I5u);
    double vcoeff = 1.0 - 1 / sqrt(I5v);

    if (I5u > 1)
    {
        ucoeff += 1.5 * strainRate * (sqrt(I5u) + 1 / sqrt(I5u) - 2);
    }
    if (I5v > 1)
    {
        vcoeff += 1.5 * strainRate * (sqrt(I5v) + 1 / sqrt(I5v) - 2);
    }




    stretch_pk1 = ucoeff * 2. * F * anisotropic_a * anisotropic_a.transpose() +
        vcoeff * 2. * F * anisotropic_b * anisotropic_b.transpose();
    return stretchS * stretch_pk1 + shearS * shear_pk1;
}

//Matrix3d computePEPF_Aniostropic3D_double(const Matrix3d &F, Vector3d direction, const double &scale,
//                                          const double &contract_length) {
//
//    direction.normalize();
//    SVDResult3D_double svdResult = QRSVD(F);
//    Matrix3d U, V, R, S, sigma;
//
//    sigma = svdResult.SIGMA;
//    V = svdResult.V;
//
//    S = V * sigma * V.transpose();
//
//    double I4 = direction.transpose() * S * direction;
//    double I5 = direction.transpose() * S.transpose() * S * direction;
//
//    if (I4 == 0) {
//        // system("pause");
//    }
//
//    double s = 0;
//    if (I4 < 0) {
//        s = -1;
//    } else if (I4 > 0) {
//        s = 1;
//    }
//
//    Matrix3d PEPF = scale * (1 - s * contract_length / sqrt(I5)) * F * direction * direction.transpose();
//    return PEPF;
//}

Matrix3d computePEPF_AniostropicRehabi3D_double(const Matrix3d &F, Vector3d direction, const double &u_anios) {
    direction.normalize();
    Matrix3d PEPF = 2 * u_anios * F * direction * direction.transpose();
    return PEPF;
}

//MatrixXd project_StabbleNHK_H_3D(const Matrix3d &F, const double &lengthRate, const double &volumRate) {
//    SVDResult3D_double svdResult = QRSVD(F);
//    Matrix3d U, sigma, V, A;
//    U = svdResult.U;
//    sigma = svdResult.SIGMA;
//    V = svdResult.V;
//
//    double I3 = sigma(0, 0) * sigma(1, 1) * sigma(2, 2);
//
//    double u = lengthRate, r = volumRate;
//    A.setZero();
//    A(0, 0) = u + r * sigma(1, 1) * sigma(2, 2) * sigma(1, 1) * sigma(2, 2);
//    A(0, 1) = sigma(2, 2) * (r * (2 * I3 - 1) - u);
//    A(0, 2) = sigma(1, 1) * (r * (2 * I3 - 1) - u);
//    A(1, 0) = sigma(2, 2) * (r * (2 * I3 - 1) - u);
//    A(1, 1) = u + r * sigma(0, 0) * sigma(2, 2) * sigma(0, 0) * sigma(2, 2);
//    A(1, 2) = sigma(0, 0) * (r * (2 * I3 - 1) - u);
//    A(2, 0) = sigma(1, 1) * (r * (2 * I3 - 1) - u);
//    A(2, 1) = sigma(0, 0) * (r * (2 * I3 - 1) - u);
//    A(2, 2) = u + r * sigma(1, 1) * sigma(0, 0) * sigma(1, 1) * sigma(0, 0);
//
//    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> eigenSolver(A);
//    //if (eigenSolver.eigenvalues()[0] >= 0.0) {
//    //    return;
//    //}
//
//    double aa = 1;
//    double bb = -(A(0, 0) + A(1, 1) + A(2, 2));
//    double cc = A(0, 0) * A(1, 1) + A(0, 0) * A(2, 2) + A(2, 2) * A(1, 1) - A(0, 1) * A(0, 1) - A(0, 2) * A(0, 2) -
//                A(2, 1) * A(2, 1);
//    double dd = A(1, 1) * A(0, 2) * A(0, 2) + A(2, 2) * A(1, 0) * A(1, 0) + A(0, 0) * A(1, 2) * A(1, 2) -
//                2 * A(1, 0) * A(1, 2) * A(0, 2) - A(0, 0) * A(1, 1) * A(2, 2);
//    double lamda[9];
//    vector<double> roots = __SolverForCubicEquation(aa, bb, cc, dd, 1e-1);
//    int rootNum = roots.size();
//    for (int i = 0; i < rootNum; i++) {
//        lamda[i] = roots[i];
//    }
//    //lamda[0] = eigenSolver.eigenvalues()[0];//A.eigenvalues()(0).real();
//    //lamda[1] = eigenSolver.eigenvalues()[1];//A.eigenvalues()(1).real();
//    //lamda[2] = eigenSolver.eigenvalues()[2];//A.eigenvalues()(2).real();
//    printf("%d  %f %f %f\n\n", rootNum, lamda[0], lamda[1], lamda[2]);
//    lamda[3] = u + sigma(2, 2) * (r * (I3 - 1) - u);
//    lamda[4] = u + sigma(0, 0) * (r * (I3 - 1) - u);
//    lamda[5] = u + sigma(1, 1) * (r * (I3 - 1) - u);
//
//    lamda[6] = u - sigma(2, 2) * (r * (I3 - 1) - u);
//    lamda[7] = u - sigma(0, 0) * (r * (I3 - 1) - u);
//    lamda[8] = u - sigma(1, 1) * (r * (I3 - 1) - u);
//    Matrix3d D0, D1, D2;
//    D0 << 1, 0, 0, 0, 0, 0, 0, 0, 0;
//    D1 << 0, 0, 0, 0, 1, 0, 0, 0, 0;
//    D2 << 0, 0, 0, 0, 0, 0, 0, 0, 1;
//
//    D0 = U * D0 * V.transpose();
//    D1 = U * D1 * V.transpose();
//    D2 = U * D2 * V.transpose();
//
//    Matrix3d Q[9];
//
//    for (int i = 0; i < 3; i++) {
//        double z0 = sigma(1, 1) * lamda[i] + sigma(0, 0) * sigma(2, 2);
//        double z1 = sigma(0, 0) * lamda[i] + sigma(1, 1) * sigma(2, 2);
//        double z2 = lamda[i] * lamda[i] - sigma(2, 2) * sigma(2, 2);
//        Q[i] = z0 * D0 + z1 * D1 + z2 * D2;
//    }
//
//    Q[3] << 0, -1, 0, 1, 0, 0, 0, 0, 0;
//    Q[4] << 0, 0, 0, 0, 0, 1, 0, -1, 0;
//    Q[5] << 0, 0, 1, 0, 0, 0, -1, 0, 0;
//    Q[6] << 0, 1, 0, 1, 0, 0, 0, 0, 0;
//    Q[7] << 0, 0, 0, 0, 0, 1, 0, 1, 0;
//    Q[8] << 0, 0, 1, 0, 0, 0, 1, 0, 0;
//
//    double ml = 1 / sqrt(2);
//    Q[3] = U * (ml * Q[3]) * V.transpose();
//    Q[4] = U * (ml * Q[4]) * V.transpose();
//    Q[5] = U * (ml * Q[5]) * V.transpose();
//    Q[6] = U * (ml * Q[6]) * V.transpose();
//    Q[7] = U * (ml * Q[7]) * V.transpose();
//    Q[8] = U * (ml * Q[8]) * V.transpose();
//    MatrixXd H = MatrixXd::Zero(9, 9);
//    //H.setZero(9, 9);
//
//    for (int i = 0; i < 9; i++) {
//        if (lamda[i] > 0) {
//            H += lamda[i] * vec_double(Q[i]) * vec_double(Q[i]).transpose();
//        }
//    }
//    return H;
//}

//MatrixXd project_StabbleNHK_2_H_3D(const Matrix3d &F, const double &lengthRate, const double &volumRate) {
//    //cout << lengthRate << "   " << volumRate << endl;
//    SVDResult3D_double svdResult = QRSVD(F);
//    Matrix3d U, sigma, V, A;
//    U = svdResult.U;
//    sigma = svdResult.SIGMA;
//    V = svdResult.V;
//
//    double I3 = sigma(0, 0) * sigma(1, 1) * sigma(2, 2);
//    double I2 = sigma(0, 0) * sigma(0, 0) + sigma(1, 1) * sigma(1, 1) + sigma(2, 2) * sigma(2, 2);
//    double g2 = sigma(0, 0) * sigma(1, 1) * sigma(0, 0) * sigma(1, 1) +
//                sigma(0, 0) * sigma(2, 2) * sigma(0, 0) * sigma(2, 2) +
//                sigma(2, 2) * sigma(1, 1) * sigma(2, 2) * sigma(1, 1);
//
//    double u = lengthRate, r = volumRate;
//    Matrix<double, 9, 9> H = MatrixXd::Zero(9, 9);
//    if (false) {
//        double n = 2 * u / ((I2 + 1) * (I2 + 1) * (r * (I3 - 1) - 3 * u / 4));
//        double p = r / (r * (I3 - 1) - 3 * u / 4);
//        double c2 = -g2 * p - I2 * n;
//        double c1 = -(1 + 2 * I3 * p) * I2 - 6 * I3 * n + (g2 * I2 - 9 * I3 * I3) * p * n;
//        double c0 = -(2 + 3 * I3 * p) * I3 + (I2 * I2 - 4 * g2) * n + 2 * I3 * p * n * (I2 * I2 - 3 * g2);
//
//        vector<double> roots = NewtonSolverForCubicEquation_snk(1, c2, c1, c0);
//
//        Matrix3d D[3];
//        double q[3];
//        Matrix3d Q[9];
//        double lamda[9];
//        double Ut = u * (1 - 1 / (I2 + 1));
//        double alpha = 1 + 3 * u / r / 4;
//
//        for (int i = 0; i < 3; i++) {
//            double alpha0 = roots[i] * (sigma(1, 1) + sigma(0, 0) * sigma(2, 2) * n + I3 * sigma(1, 1) * p) +
//                            sigma(0, 0) * sigma(2, 2) + sigma(1, 1) *
//                                                        (sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1) +
//                                                         sigma(2, 2) * sigma(2, 2)) * n +
//                            I3 * sigma(0, 0) * sigma(2, 2) * p +
//                            sigma(0, 0) * (sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1)) * sigma(2, 2) *
//                            (sigma(1, 1) * sigma(1, 1) - sigma(2, 2) * sigma(2, 2)) * p * n;
//
//            double alpha1 = roots[i] * (sigma(0, 0) + sigma(1, 1) * sigma(2, 2) * n + I3 * sigma(0, 0) * p) +
//                            sigma(1, 1) * sigma(2, 2) - sigma(0, 0) *
//                                                        (sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1) -
//                                                         sigma(2, 2) * sigma(2, 2)) * n +
//                            I3 * sigma(1, 1) * sigma(2, 2) * p -
//                            sigma(1, 1) * (sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1)) * sigma(2, 2) *
//                            (sigma(0, 0) * sigma(0, 0) - sigma(2, 2) * sigma(2, 2)) * p * n;
//
//            double alpha2 = roots[i] * roots[i] - roots[i] * (sigma(0, 0) * sigma(0, 0) + sigma(1, 1) * sigma(1, 1)) *
//                                                  (n + sigma(2, 2) * sigma(2, 2) * p) -
//                            sigma(2, 2) * sigma(2, 2) - 2 * I3 * n - 2 * I3 * sigma(2, 2) * sigma(2, 2) * p +
//                            ((sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1)) * sigma(2, 2)) *
//                            ((sigma(0, 0) * sigma(0, 0) - sigma(1, 1) * sigma(1, 1)) * sigma(2, 2)) * p * n;
//
//            double CheckZeroV = alpha0 * alpha0 + alpha1 * alpha1 + alpha2 * alpha2;
//
//            if (CheckZeroV == 0) {
//                lamda[i] = 0;
//                continue;
//            }
//
//            q[i] = 1 / sqrt(CheckZeroV);
//            D[i] << alpha0, 0, 0, 0, alpha1, 0, 0, 0, alpha2;
//            Q[i] = q[i] * U * D[i] * V.transpose();
//            lamda[i] = r * (I3 - alpha) * roots[i] + Ut;
//        }
//
//        lamda[3] = Ut + sigma(2, 2) * r * (I3 - alpha);
//        lamda[4] = Ut + sigma(0, 0) * r * (I3 - alpha);
//        lamda[5] = Ut + sigma(1, 1) * r * (I3 - alpha);
//
//        lamda[6] = Ut - sigma(2, 2) * r * (I3 - alpha);
//        lamda[7] = Ut - sigma(0, 0) * r * (I3 - alpha);
//        lamda[8] = Ut - sigma(1, 1) * r * (I3 - alpha);
//
//        Q[3] << 0, -1, 0, 1, 0, 0, 0, 0, 0;
//        Q[4] << 0, 0, 0, 0, 0, 1, 0, -1, 0;
//        Q[5] << 0, 0, 1, 0, 0, 0, -1, 0, 0;
//        Q[6] << 0, 1, 0, 1, 0, 0, 0, 0, 0;
//        Q[7] << 0, 0, 0, 0, 0, 1, 0, 1, 0;
//        Q[8] << 0, 0, 1, 0, 0, 0, 1, 0, 0;
//
//        double ml = 1 / sqrt(2);
//        Q[3] = ml * U * Q[3] * V.transpose();
//        Q[4] = ml * U * Q[4] * V.transpose();
//        Q[5] = ml * U * Q[5] * V.transpose();
//        Q[6] = ml * U * Q[6] * V.transpose();
//        Q[7] = ml * U * Q[7] * V.transpose();
//        Q[8] = ml * U * Q[8] * V.transpose();
//
//        for (int i = 0; i < 9; i++) {
//            if (lamda[i] > 0) {
//                H += lamda[i] * vec_double(Q[i]) * vec_double(Q[i]).transpose();
//            }
//        }
//    } else {
//        double mu = u, lambda = r;
//        double Ic = I2;
//        Eigen::Matrix<double, 9, 9> H1 = 2 * Eigen::Matrix<double, 9, 9>::Identity();
//        Eigen::Matrix<double, 9, 1> g1;
//        g1.block<3, 1>(0, 0) = 2 * F.col(0);
//        g1.block<3, 1>(3, 0) = 2 * F.col(1);
//        g1.block<3, 1>(6, 0) = 2 * F.col(2);
//        Eigen::Matrix<double, 9, 1> gJ;
//        gJ.block<3, 1>(0, 0) = F.col(1).cross(F.col(2));
//        gJ.block<3, 1>(3, 0) = F.col(2).cross(F.col(0));
//        gJ.block<3, 1>(6, 0) = F.col(0).cross(F.col(1));
//        Eigen::Matrix3d f0hat;
//        f0hat <<
//              0, -F(2, 0), F(1, 0),
//                F(2, 0), 0, -F(0, 0),
//                -F(1, 0), F(0, 0), 0;
//        Eigen::Matrix3d f1hat;
//        f1hat <<
//              0, -F(2, 1), F(1, 1),
//                F(2, 1), 0, -F(0, 1),
//                -F(1, 1), F(0, 1), 0;
//        Eigen::Matrix3d f2hat;
//        f2hat <<
//              0, -F(2, 2), F(1, 2),
//                F(2, 2), 0, -F(0, 2),
//                -F(1, 2), F(0, 2), 0;
//        Eigen::Matrix<double, 9, 9> HJ;
//        HJ.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
//        HJ.block<3, 3>(0, 3) = -f2hat;
//        HJ.block<3, 3>(0, 6) = f1hat;
//        HJ.block<3, 3>(3, 0) = f2hat;
//        HJ.block<3, 3>(3, 3) = Eigen::Matrix3d::Zero();
//        HJ.block<3, 3>(3, 6) = -f0hat;
//        HJ.block<3, 3>(6, 0) = -f1hat;
//        HJ.block<3, 3>(6, 3) = f0hat;
//        HJ.block<3, 3>(6, 6) = Eigen::Matrix3d::Zero();
//        double J = I3;
//        H = (Ic * mu) / (2 * (Ic + 1)) * H1 + lambda * (J - 1 - (3 * mu) / (4.0 * lambda)) * HJ +
//            (mu / (2 * (Ic + 1) * (Ic + 1))) * g1 * g1.transpose() + lambda * gJ * gJ.transpose();
//        IglUtils::makePD<double, 9>(H);
//    }
//    return H;
//}

Eigen::Matrix<double, 6, 6> project_baraffwitkint_H_3D(const Matrix<double, 3, 2>& F,
    const Vector2d& anisotropic_a,
    const Vector2d& anisotropic_b,
    double stretchS, double shearS, double strainRate) {
    Eigen::Matrix<double, 6, 6> final_H = Eigen::Matrix<double, 6, 6>::Zero();
    {
        Eigen::Matrix<double, 6, 6> H;
        H.setZero();
        double I5u = (F * anisotropic_a).transpose() * F * anisotropic_a;
        double I5v = (F * anisotropic_b).transpose() * F * anisotropic_b;
        double invSqrtI5u = 1.0 / sqrt(I5u);
        double invSqrtI5v = 1.0 / sqrt(I5v);

        double sqrtI5u = sqrt(I5u);
        double sqrtI5v = sqrt(I5v);

        if (sqrtI5u > 1)
            H(0, 0) = H(1, 1) = H(2, 2) = 2 * (((sqrtI5u - 1) * (3 * sqrtI5u * strainRate - 3 * strainRate + 2)) / (2 * sqrtI5u));
        if (sqrtI5v > 1)
            H(3, 3) = H(4, 4) = H(5, 5) = 2 * (((sqrtI5v - 1) * (3 * sqrtI5v * strainRate - 3 * strainRate + 2)) / (2 * sqrtI5v));
        auto fu = F.col(0).normalized();
        auto fv = F.col(1).normalized();

        double uCoeff = (sqrtI5u > 1.0) ? (3 * I5u * strainRate - 3 * strainRate + 2)
            / (sqrt(I5u)) :
            2.0;
        double vCoeff = (sqrtI5v > 1.0) ? (3 * I5v * strainRate - 3 * strainRate + 2)
            / (sqrt(I5v)) :
            2.0;


        H.block<3, 3>(0, 0) += uCoeff * (fu * fu.transpose());
        H.block<3, 3>(3, 3) += vCoeff * (fv * fv.transpose());

        final_H += stretchS * H;
    }
    {
        Eigen::Matrix<double, 6, 6> H_shear;
        H_shear.setZero();
        Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
        H(3, 0) = H(4, 1) =
            H(5, 2) = H(0, 3) =
            H(1, 4) = H(2, 5) = 1.0;
        double I6 = anisotropic_a.transpose() * F.transpose() * F * anisotropic_b;
        double signI6 = (I6 >= 0) ? 1.0 : -1.0;
        auto g = F * (anisotropic_a * anisotropic_b.transpose() + anisotropic_b * anisotropic_a.transpose());
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
        final_H += shearS * H_shear;
    }
    return final_H;

}

//MatrixXd
//project_ANIOSI5_H_3D(const Matrix3d &F, Vector3d direction, const double &scale, const double &contract_length) {
//    direction.normalize();
//    SVDResult3D_double svdResult = QRSVD(F);
//    Matrix3d U, sigma, V, S;
//    U = svdResult.U;
//    sigma = svdResult.SIGMA;
//    V = svdResult.V;
//
//    S = V * sigma * V.transpose();
//    double I4 = direction.transpose() * S * direction;
//    double I5 = direction.transpose() * S.transpose() * S * direction;
//
//    if (abs(I5) < 1e-15) return MatrixXd::Zero(9, 9);
//
//    double s = 0;
//    if (I4 < 0) {
//        s = -1;
//    } else if (I4 > 0) {
//        s = 1;
//    }
//
//    double lamda0 = scale;
//    double lamda1 = scale * (1 - s * contract_length / sqrt(I5));
//    double lamda2 = lamda1;
//    //double lamda2 = lamda1;
//    Matrix3d Q0, Q1, Q2, A;
//    A = direction * direction.transpose();
//    Q0 = (1 / sqrt(I5)) * F * A;
//
//    Matrix3d Tx, Ty, Tz;
//
//    Tx << 0, 0, 0, 0, 0, 1, 0, -1, 0;
//    Ty << 0, 0, -1, 0, 0, 0, 1, 0, 0;
//    Tz << 0, 1, 0, -1, 0, 0, 0, 0, 0;
//    Vector3d directionM = V.transpose() * direction;
//
//    Tx *= 1.f / sqrt(2.f);
//    Ty *= 1.f / sqrt(2.f);
//    Tz *= 1.f / sqrt(2.f);
//
//    Q1 = U * Tx * sigma * V.transpose() * A;
//    Q2 = (sigma(1, 1) * directionM[1]) * U * Tz * sigma * V.transpose() * A -
//         (sigma(2, 2) * directionM[2]) * U * Ty * sigma * V.transpose() * A;
//
//    MatrixXd H = lamda0 * vec_double(Q0) * vec_double(Q0).transpose();
//    if (lamda1 > 0) {
//        H += lamda1 * vec_double(Q1) * vec_double(Q1).transpose();
//        H += lamda2 * vec_double(Q2) * vec_double(Q2).transpose();
//    }
//
//    return H;
//}

//MatrixXd project_ANIOSI5_Rehabi_H_3D(const Matrix3d &F, Vector3d direction, const double &u_anios) {
//    direction.normalize();
//    SVDResult3D_double svdResult = QRSVD(F);
//    Matrix3d U, sigma, V, S;
//    U = svdResult.U;
//    sigma = svdResult.SIGMA;
//    V = svdResult.V;
//
//    S = V * sigma * V.transpose();
//    //double I4 = direction.transpose() * S * direction;
//
//
//    //if (abs(I5) < 1e-15) return MatrixXd::Zero(9, 9);
//
//    //double s = 0;
//    //if (I4 < 0) {
//    //    s = -1;
//    //}
//    //else if (I4 > 0) {
//    //    s = 1;
//    //}
//
//    double lamda0 = 2 * u_anios;
//    double lamda1 = lamda0;
//    double lamda2 = lamda1;
//    //double lamda2 = lamda1;
//    Matrix3d Q0, Q1, Q2, A;
//
//    /*double I5 = direction.transpose() * S.transpose() * S * direction;
//    A = direction * direction.transpose();
//    Q0 = (1 / sqrt(I5)) * F * A;
//
//    Matrix3d Tx, Ty, Tz;
//
//    Tx << 0, 0, 0, 0, 0, 1, 0, -1, 0;
//    Ty << 0, 0, -1, 0, 0, 0, 1, 0, 0;
//    Tz << 0, 1, 0, -1, 0, 0, 0, 0, 0;
//    Vector3d directionM = V.transpose() * direction;
//
//    Tx *= 1.f / sqrt(2.f);
//    Ty *= 1.f / sqrt(2.f);
//    Tz *= 1.f / sqrt(2.f);
//
//    Q1 = U * Tx * sigma * V.transpose() * A;
//    Q2 = (sigma(1, 1) * directionM[1]) * U * Tz * sigma * V.transpose() * A - (sigma(2, 2) * directionM[2]) * U * Ty * sigma * V.transpose() * A;*/
//
//    Q0 << direction[0], direction[1], direction[2], 0, 0, 0, 0, 0, 0;
//    Q1 << 0, 0, 0, direction[0], direction[1], direction[2], 0, 0, 0;
//    Q2 << 0, 0, 0, 0, 0, 0, direction[0], direction[1], direction[2];
//
//
//    MatrixXd H = lamda0 * vec_double(Q0) * vec_double(Q0).transpose();
//    H += lamda1 * vec_double(Q1) * vec_double(Q1).transpose();
//    H += lamda2 * vec_double(Q2) * vec_double(Q2).transpose();
//
//    return H;
//}

}
