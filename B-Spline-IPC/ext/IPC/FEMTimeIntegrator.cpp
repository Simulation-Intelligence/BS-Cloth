#include "FEMTimeIntegrator.h"
#include "fem_parameters.h"
//#include "fem3D_Device.cuh"

namespace IPC
{

    ImplicitFEMIntegrator::ImplicitFEMIntegrator(model_tet* tetra_mesh3d, unsigned int m_sceneType) {
        meshTetes = tetra_mesh3d;
        sceneType = m_sceneType;
    }



    int ImplicitFEMIntegrator::integrate(int& stepId, int& total_cg_iterations, int& total_newton_iterations, SpatialHash& sh, Ground& gd) {

        mesh3D& mesh = meshTetes->mesh3Ds[0];

        return IPC_Solver(stepId, mesh, sh, gd);
    }
}





