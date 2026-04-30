#include "fem_parameters.h"
#include "Simulator.h"
#include "iostream"
#include "Eigen/Eigen"
#include <fstream>

namespace IPC
{

    using namespace Eigen;
    using namespace FEM;

    void buildSpecialPoints(mesh3D& mesh) {
        Vector3d tp[6];
        tp[0] = Vector3d(-0.26554, 0.67295, -0.23324);
        tp[1] = Vector3d(-0.05724, 0.69143, -0.10306);
        tp[2] = Vector3d(0.04884, 0.47191, 0.04952);
        tp[3] = Vector3d(-0.26145, -0.12829, -0.23713);
        tp[4] = Vector3d(-0.06139, -0.10531, -0.09838);
        tp[5] = Vector3d(0.015541, -0.32489, 0.02275);
        mesh.spectialPontsArray.resize(6);
        for (int i = 0; i < mesh.spectialPontsArray.size(); i++) {
            double mindist = 1e32;
            for (int j = 0; j < mesh.surfVerts.size(); j++) {
                double dist = (mesh.vertices[mesh.surfVerts[j]] - tp[i]).norm();
                //printf("dist: %f\n", dist);
                if (dist < mindist) {
                    mindist = dist;
                    mesh.spectialPontsArray[i] = mesh.surfVerts[j];
                    //printf("special Point: %d\n", mesh.surfVerts[j]);
                }
            }
            //printf("special Point: %d\n", mesh.spectialPontsArray[i]);
            //cout<<mesh.vertexes[mesh.spectialPontsArray[i]]<<endl;
        }


    }

    void DefaultSettings(mesh3D& mesh3d) {

        mesh3d.density = 1e3;
        mesh3d.cloth_density = 1e2;
        mesh3d.clothThickness = 1e-3;

        mesh3d.use_barrier = 1;
        mesh3d.IPC_dt = 1 / 240.; // 240 fps

        mesh3d.drag_coeff = 0.00;
        mesh3d.Hhat = 9e-8;

        mesh3d.Fhat = 1e-6;
        mesh3d.Kappa = 0;
        mesh3d.dTol = 1e-18;

        mesh3d.YoungModulus = 1e4;
        mesh3d.PoissonRate = 0.49;
        mesh3d.friction = 0.5;

        mesh3d.clothYoungModulus = 1e4;
        mesh3d.bendYoungModulus = 1e6;
        mesh3d.strainRate = 1e2;

        mesh3d.Newton_Solver_Threshold = 1e-2;

        mesh3d.lengthRateLame = mesh3d.YoungModulus / (2 * (1 + mesh3d.PoissonRate));
        mesh3d.volumeRateLame = mesh3d.YoungModulus * mesh3d.PoissonRate / ((1 + mesh3d.PoissonRate) * (1 - 2 * mesh3d.PoissonRate));
        mesh3d.lengthRate = 4 * mesh3d.lengthRateLame / 3;
        mesh3d.volumeRate = mesh3d.volumeRateLame + 5 * mesh3d.lengthRateLame / 6;
        mesh3d.stretchStiffness = mesh3d.clothYoungModulus / (2 * (1 + mesh3d.PoissonRate));
        //mesh3d.shearStiffness = mesh3d.stretchStiffness * 1;
        mesh3d.shearStiffness = mesh3d.stretchStiffness * 1;
        mesh3d.bendingStiffness = mesh3d.bendYoungModulus * pow(mesh3d.clothThickness, 3) / (24 * (1 - mesh3d.PoissonRate * mesh3d.PoissonRate));

    }


    void LoadSettings(mesh3D& mesh3d) {
        bool successfulRead = false;

        //read file
        std::ifstream infile;


#ifdef _WIN32
        string DEFAULT_CONFIG_FILE = "scene/parameterSetting.txt";
#else
        string DEFAULT_CONFIG_FILE = "../CPU IPC/scene/parameterSetting.txt";
#endif


        infile.open(DEFAULT_CONFIG_FILE, std::ifstream::in);
        if (successfulRead = infile.is_open())
        {
            char ignoreToken[256];
            mesh3d.Fhat = 1e-6;
            mesh3d.Kappa = 0;
            mesh3d.dTol = 1e-18;

            // global settings:
            infile >> ignoreToken >> mesh3d.density;
            infile >> ignoreToken >> mesh3d.YoungModulus;
            infile >> ignoreToken >> mesh3d.PoissonRate;
            infile >> ignoreToken >> mesh3d.friction;
            infile >> ignoreToken >> mesh3d.clothThickness;
            infile >> ignoreToken >> mesh3d.clothYoungModulus;
            infile >> ignoreToken >> mesh3d.bendYoungModulus;
            infile >> ignoreToken >> mesh3d.cloth_density;
            infile >> ignoreToken >> mesh3d.strainRate;
            infile >> ignoreToken >> mesh3d.use_barrier;
            infile >> ignoreToken >> mesh3d.drag_coeff;
            infile >> ignoreToken >> mesh3d.IPC_dt;
            infile >> ignoreToken >> mesh3d.Newton_Solver_Threshold;
            infile >> ignoreToken >> mesh3d.Hhat;

            mesh3d.Hhat *= mesh3d.Hhat;
            mesh3d.lengthRateLame = mesh3d.YoungModulus / (2 * (1 + mesh3d.PoissonRate));
            mesh3d.volumeRateLame = mesh3d.YoungModulus * mesh3d.PoissonRate / ((1 + mesh3d.PoissonRate) * (1 - 2 * mesh3d.PoissonRate));
            mesh3d.lengthRate = 4 * mesh3d.lengthRateLame / 3;
            mesh3d.volumeRate = mesh3d.volumeRateLame + 5 * mesh3d.lengthRateLame / 6;
            mesh3d.stretchStiffness = mesh3d.clothYoungModulus / (2 * (1 + mesh3d.PoissonRate));
            //mesh3d.shearStiffness = mesh3d.stretchStiffness * 0.0003;
            mesh3d.shearStiffness = 1000.;
            mesh3d.bendingStiffness = mesh3d.clothYoungModulus * pow(mesh3d.clothThickness, 3) / (24 * (1 - mesh3d.PoissonRate * mesh3d.PoissonRate));

            infile.close();
        }

        if (!successfulRead)
        {
            std::cerr << "Waning: failed loading settings, set to defaults." << std::endl;
            DefaultSettings(mesh3d);
        }
    }


    bool FEMSimulator::buildModels(unsigned int buildType, unsigned int sceneType) {
        mesh3D mesh3d;
        LoadSettings(mesh3d);

        mesh3d.maxCorner = Vector3d(-1e32, -1e32, -1e32);
        mesh3d.minCorner = Vector3d(1e32, 1e32, 1e32);

        mesh3d.objMaxConer = Vector3d(0, 0, 0);
        mesh3d.objMinConer = Vector3d(0, 0, 0);

        if (true)
        {
            //mesh3d.load_triangleMesh("../CPU IPC/triangleMesh/tricloth.obj", 0.3, Vector3d(0, 0, 0));
            //mesh3d.load_triangleMesh("../CPU IPC/triangleMesh/CMU/plane100.obj", 0.3, Vector3d(0, 0, 0));

            mesh3d.load_triangleMesh("../CPU IPC/bsipc-test/mesh-49.obj", 1.0, Vector3d(0, 0, 0));

            //Matrix3d rotate, rotatey;
            //float angleX = FEM::PI / 2, angleY = FEM::PI / 4, angleZ = FEM::PI / 2;
            //rotate << 1, 0, 0, 0, cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);
            //rotatey << cos(angleY), -sin(angleY), 0, sin(angleY), cos(angleY), 0, 0, 0, 1;
            //for (int j = 0; j < mesh3d.vertexNum; j++) {
            //    mesh3d.vertexes[j] = (rotate * mesh3d.vertexes[j]);

            //    //mesh3d.vertexes[j] = (rotatey * mesh3d.vertexes[j]);
            //}
            mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0, 0));
            //mesh3d.load_triangleMesh("../CPU IPC/triangleMesh/tricloth.obj", 0.7, Vector3d(0, 0, 0));
            double xmin = 1e32, ymin = 1e32, zmin = 1e32;
            double xmax = -1e32, ymax = -1e32, zmax = -1e32;
            for (int j = 0; j < mesh3d.vertexNum; j++) {

                Vector3d pos = mesh3d.vertices[j];
                if (xmin > pos[0]) xmin = pos[0];
                if (ymin > pos[1]) ymin = pos[1];
                if (zmin > pos[2]) zmin = pos[2];
                if (xmax < pos[0]) xmax = pos[0];
                if (ymax < pos[1]) ymax = pos[1];
                if (zmax < pos[2]) zmax = pos[2];


            }
            mesh3d.maxCorner = Vector3d(xmax, ymax, zmax);
            mesh3d.minCorner = Vector3d(xmin, ymin, zmin);

            mesh3d.boundaryTypes[0] = 1;
            mesh3d.boundaryTypes[49] = 1;
            mesh3d.Constraints[0].setZero();
            mesh3d.Constraints[49].setZero();

        }

        //mesh3d.load_tetrahedraMesh_IPC_TetMesh("../CPU IPC/tetrahedraMesh/ipcmesh/cube.msh", 0.5, Vector3d(0.25, -0.6, 0.25));
        //mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny.msh", 0.2, Vector3d(0, -0.5, 0));
        //mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0.15, 0));
        //mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny.msh", 0.5, Vector3d(0, 0.0, 0));
        //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", 0.5, Vector3d(0.25, 0.6, 0.25));
        //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", 0.5, Vector3d(0.25, 0.8, 0.25));
        //mesh3d.load_tetrahedraMesh("tetrahedraMesh/cube15.msh", 0.5, Vector3d(-0.25, -0.95, -0.25));
        //mesh3d.load_tetrahedraMesh("tetrahedraMesh/ipcmesh/sqballTet_.msh", 1, Vector3d(0, 0.75, 0));
        //mesh3d.load_tetrahedraMesh("tetrahedraMesh/ipcmesh/sqballTet_.msh", 1, Vector3d(0, -0.3, 0));
        //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/sqballTet_.msh", 1, Vector3d(0, -0.75, 0));


        mesh3d.vertexNum = mesh3d.vertices.size();
        mesh3d.tetrahedralNum = mesh3d.tetrahedrals.size();
        mesh3d.triangleNum = mesh3d.triangles.size();

        //mesh3d.triangleNum = 0;
        printf("tets num: %d\n", mesh3d.tetrahedralNum);
        //printf("triangles num: %d\n", mesh3d.triangleNum);
        printf("verts num: %d\n", mesh3d.vertexNum);

        initMesh3D(mesh3d, 1, 0.2);
        //mesh3d.triangleNum = 0;

        //for (int i = 0; i < mesh3d.tetrahedraNum; i++) {
        //    Vector4i tet = Vector4i(mesh3d.tetrahedras[i][0], mesh3d.tetrahedras[i][1], mesh3d.tetrahedras[i][2], mesh3d.tetrahedras[i][3]);
        //    mesh3d.tetrahedrasV.push_back(tet);
        //}

        //mesh3d.restSNKE = getObjRestEnergy_StableNHK2_3D(mesh3d.vertexes, mesh3d, FEM::lengthRate, FEM::volumeRate) * FEM::IPC_dt * FEM::IPC_dt;

        mesh3d.v_rest = mesh3d.vertices;
        mesh3d.V_prev = mesh3d.vertices;

        ComputeXTilde(mesh3d);


        mesh3d.bboxDiagSize2 = (mesh3d.maxCorner - mesh3d.minCorner).squaredNorm();
        mesh3d.Hhat *= mesh3d.bboxDiagSize2;// (mesh3d.objMaxConer - mesh3d.objMinConer).squaredNorm();
        mesh3d.Fhat *= mesh3d.bboxDiagSize2;;
        mesh3d.dTol *= mesh3d.bboxDiagSize2;


        tetrahedra_meshes.mesh3Ds.push_back(mesh3d);
        tetrahedra_meshes.calculate_surface();
        //tetrahedra_meshes.mesh3Ds[0].triangles.clear();
        //buildSpecialPoints(tetrahedra_meshes.mesh3Ds[0]);

        double length = 0;
        for (const auto& edg : (tetrahedra_meshes.mesh3Ds[0]).surfEdges) {
            length += (tetrahedra_meshes.mesh3Ds[0].vertices[edg.first] -
                tetrahedra_meshes.mesh3Ds[0].vertices[edg.second]).norm();
        }
        tetrahedra_meshes.mesh3Ds[0].averageEdgeLength = length / (tetrahedra_meshes.mesh3Ds[0].surfEdges.size() * 3);

        printf("triangles num: %d\n", tetrahedra_meshes.mesh3Ds[0].surface.size());
        buildIntegrator(buildType, sceneType);

        if (!tetrahedra_meshes.mesh3Ds[0].load_tetTempData()) {
            printf("no temp data\n");
        }
        else {
            printf("load temp data\n");
        }

        //if (tetrahedra_meshes.mesh3Ds[0].use_barrier) {

        BuildCollisionSets();
        //}
        //printf("finish build cp\n");

        return true;
    }


    void FEMSimulator::buildIntegrator(const int integratorType, unsigned int sceneType) {
        integrator = std::make_unique<ImplicitFEMIntegrator>(&tetrahedra_meshes, sceneType);
    }

    void FEMSimulator::BuildCollisionSets() {
        if (tetrahedra_meshes.mesh3Ds[0].use_barrier) {
            sh.build(tetrahedra_meshes.mesh3Ds[0], tetrahedra_meshes.mesh3Ds[0].averageEdgeLength);
            sh.calculateActiveSet(tetrahedra_meshes.mesh3Ds[0]);
        }

        gd.calculateActiveSet(tetrahedra_meshes.mesh3Ds[0]);
    }

    int FEMSimulator::simulateStick(int& stepId) {
        int cg_loops = 0;
        int newTon_loops = 0;
        int k = integrator->integrate(stepId, cg_loops, newTon_loops, sh, gd);
        //buildCollisionSets();
        return k;
    }

}