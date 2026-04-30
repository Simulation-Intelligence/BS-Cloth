#pragma once
#ifndef _FEM_TIME_INTEGRATOR_CUH_
#define _FEM_TIME_INTEGRATOR_CUH_

#include "mesh.h"
//#include "FEMMeshes.cuh"
#include <memory>
//#include "Device_PCG.cuh"
#include "collisionUtil.h"
#include "IPC_FUNC.h"

namespace IPC
{

	class FEMIntegrator {
	public:
		FEMIntegrator() {};
		FEMIntegrator(int vertexNUm, int tetrahedralNum);
		~FEMIntegrator() {};
		virtual int integrate(int& stepId, int& total_cg_iterations, int& total_newton_iterations, SpatialHash& sh, Ground& gd) = 0;
		void solveCollision();
	public:
		int vertex_Num;
		int tetrahedra_Num;
		int modified_Num;
		int sceneType;
		int loopCount;
	};



	class ImplicitFEMIntegrator : public FEMIntegrator {
		friend class FEMSimulator;
	public:
		ImplicitFEMIntegrator(model_tet* tetra_mesh3d, unsigned int sceneType);
		~ImplicitFEMIntegrator() {}
		virtual int integrate(int& stepId, int& total_cg_iterations, int& total_newton_iterations, SpatialHash& sh, Ground& gd);
	private:

		//void Projected_Newton3D(mesh3D& mesh, double& mfsum, int& total_cg_iterations, int& total_newton_iterations);
		model_tet* meshTetes;
		double mfsum;
		int total_cg_iterations;
		int total_newton_iterations;


	};

}

#endif //_FEM_TIME_INTEGRATOR_CUH_