#pragma once
#ifndef FEM_MESH_H
#define FEM_MESH_H

#include <cstdint>
#include <vector>
#include <sstream>
#include <fstream>
#include <unordered_set>
#include <utility>

#include "Eigen/Eigen"
#include "mIPC.h"

namespace IPC
{

	using namespace std;
	using namespace Eigen;
	struct mesh_circle {
		double vertex[9][2] = { {0,0},
			{1,0},
		{0.7071067811865476, 0.7071067811865476},
		{0, 1},
		{-0.7071067811865476, 0.7071067811865476},
		{-1, 0},
		{-0.7071067811865476, -0.7071067811865476},
		{0, -1},
		{0.7071067811865476, -0.7071067811865476} };

		//double force[9][2] = {0};
		//double velocity[9][2] = {0};

		int index[8][3] = { {1, 2, 0},
		{2, 3, 0},
		{3, 4, 0},
		{4, 5, 0},
		{5, 6, 0},
		{6, 7, 0},
		{7, 8, 0},
		{8, 1, 0} };
		//Matrix2d DM_tetrahedra_inverse[8];
		int vertexNum = 9;
		int triangleNum = 8;
	};

	struct mesh_rectangle {
		double vertex[15][2] = { {0,1},
			{1,1},
			{1, 2},
			{0, 2},
			{-1, 2},
			{-1, 1},
			{-1, 0},
			{0, 0},
			{1, 0},
			{-1,-1},
			{0,-1},
			{1,-1},
			{-1,-2},
			{0, -2},
			{1,-2} };

		//double force[15][2] = { 0 };
		//double velocity[15][2] = { 0 };

		int index[16][3] = { {1, 2, 0},
			{2, 3, 0},
			{3, 4, 0},
			{4, 5, 0},
			{5, 6, 0},
			{6, 7, 0},
			{7, 8, 0},
			{8, 1, 0},
			{6, 7, 10},
			{7, 8, 10},
			{8, 11, 10},
			{11, 14, 10},
			{13, 14, 10},
			{12, 13, 10},
			{9, 12, 10},
			{6, 9, 10} };
		//Matrix2d DM_tetrahedra_inverse[16];
		int vertexNum = 15;
		int triangleNum = 16;
	};

	struct mesh_cuboid {
		double vertex[8][3] = {
			{-1,1,1},
			{-1,1,-1},
			{1,1,-1},
			{1,1,1},
			{-1,-1,1},
			{-1,-1,-1},
			{1,-1,-1},
			{1,-1,1} };

		//double force[8][3] = { 0 };
		//double velocity[8][3] = { 0 };

		uint64_t index[5][4] = { {0, 4, 7, 5},
			{2, 3, 7, 0},
			{2, 6, 7, 5},
			{0, 1, 2, 5},
			{0, 2, 7, 5} };
		//Matrix3d DM_tetrahedra_inverse[5];
		int vertexNum = 8;
		int tetrahedralNum = 5;
	};


	class mesh2D {
	public:
		vector<double> areas;
		vector<double> masses;
		vector<Vector2d> vertexes;
		vector<Vector3i> triangles;
		vector<Vector2d> forces;
		vector<Vector2d> velocities;
		vector<Matrix2d> DM_triangle_inverse;
		int vertexNum;
		int triangleNum;
		void InitMesh(int type, double scale);
	};

	struct SurfaceInfo
	{
		vector<size_t> surfVerts;
		vector<size_t> surfEdges;
		vector<size_t> surfaces;
	};

	class mesh3D {
	public:
		double density;
		double lengthRateLame;
		double volumeRateLame;
		double lengthRate;
		double volumeRate;
		double friction;

		double cloth_density;
		double YoungModulus;
		double PoissonRate;

		double clothThickness;
		double clothYoungModulus;
		double bendYoungModulus;
		double stretchStiffness;
		double shearStiffness;
		double strainRate;
		double bendingStiffness;
		double maxVolum;
		vector<Matrix3d> Constraints;
		//vector<bool> isDelete;
		//vector<int> idmap;
		//vector<int> idinverseMap;
		vector<double> volum;
		vector<double> areas;
		vector<double> masses;
		double meanMass;
		double restSNKE;
		double averageEdgeLength;
		double Hhat;
		double Fhat;
		double drag_coeff;
		double IPC_dt;
		double Kappa;
		double bboxDiagSize2;
		double dTol;
		double Newton_Solver_Threshold;
		bool use_barrier;

		vector<Vector3d> vertices;
		vector<Vector3d> v_rest;
		vector<Vector4i> tetrahedrals;
		vector<Vector3i> triangles;
		vector<Vector3d> velocities;
		vector<Matrix3d> DM_tetrahedra_inverse;
		vector<Matrix2d> DM_triangle_inverse;
		vector<Vector4i> surface;
		vector<uint64_t> surfVerts;
		vector<pair<uint64_t, uint64_t>> surfEdges;

		vector<double> Self_lambda_lastH;
		vector<double> Environment_lambda_lastH;
		vector<Eigen::Vector2d> MMDistCoord;
		vector<Eigen::Matrix<double, 3, 2>> MMTanBasis;
		vector<int> Environment_ActiveSet, Environment_activeSet_lastH;
		vector<MMCVID> Self_ActiveSet, Self_activeSet_lastH;
		vector<MMCVID> Self_EE_ActiveSet;
		vector<pair<int, int>> Self_EEeIe_ActiveSet;
		vector<pair<int, int>> Self_CCD_ActiveSet;
		vector<int> closeConstraintID;
		::std::vector<MMCVID> closeMConstraintID;
		vector<double> closeConstraintVal;
		vector<double> closeMConstraintVal;

		Vector3d minCorner, maxCorner;
		Vector3d objMinConer, objMaxConer;
		vector<int> boundaryTypes;
		vector<int> spectialPontsArray;
		//IPC
		vector<Vector3d> xTilta, dx_Elastic, acceleration;
		vector<Vector3d> V_prev;

		vector<Eigen::VectorXd> EKF;
		bool characterAnimation;


		vector<Vector2i> tri_edges_adj_points;
		vector<Vector2i> tri_edges;

		vector<SurfaceInfo> surfBackLinks;
        vector<size_t> patchSurfaceSep;				// the starting index of each patch in surfBackLinks
        vector<size_t> patchVertexSep;              // the starting index of each patch in surfVerts
        vector<size_t> patchEdgeSep;                // the starting index of each patch in surfEdges

		int D12x12Num;
		int D9x9Num;
		int D6x6Num;
		int D3x3Num;
		int vertexNum;
		int tetrahedralNum;
		int triangleNum;
		void InitMesh(int type, double scale);
		bool load_tetrahedraMesh(const ::std::string& filename, double scale, Vector3d offset);
		bool load_triangleMesh(const ::std::string& filename, double scale, Vector3d position_offset, int type = 0);
		bool load_tetrahedraMesh_IPC_TetMesh(const ::std::string& filename, double scale, Vector3d position_offset);
		void load_test(double scale, int num = 1);
		void getSurface();
		bool output_tetrahedraMesh(const ::std::string& filename);
		bool output_tetTempData();
		bool load_tetTempData();

		vector<Eigen::VectorXd> get_dXn1_dXn();
		vector<Eigen::VectorXd> get_dXn1_dVn();
		vector<Eigen::VectorXd> get_dVn1_dXn();
		vector<Eigen::VectorXd> get_dVn1_dVn();

        bool filterPointTrig(const unsigned int surfVertIdx, const unsigned int surfTrigIdx) const;					// returns true if the point-triangle pair should be filtered (not included in contact)
        bool filterEdgeEdge(const unsigned int surfEdgeIdxI, const unsigned int surfEdgeIdxJ) const;				// returns true if the edge-edge pair should be filtered (not included in contact)
		bool filterEdgeTrig(const unsigned int surfEdgeIdx, const unsigned int surfTrigIdx) const;					// returns true if the edge-triangle pair should be filtered (not included in contact)

		// compress the keys before inserting into the sets
        unordered_set<uint64_t> pointTrigFilterSet;
        unordered_set<uint64_t> edgeEdgeFilterSet;
        unordered_set<uint64_t> edgeTrigFilterSet;
	};

	inline uint64_t compressPair(const unsigned int a, const unsigned int b)
	{
		return (static_cast<uint64_t>(a) << 32) | b;
	}

	inline array<unsigned int, 2> decompressPair(const uint64_t key)
	{
		return { static_cast<unsigned int>(key >> 32), static_cast<unsigned int>(key & 0xFFFFFFFF) };
    }


	inline void saveSurfaceMesh(const string& path, const mesh3D& meshTemp) {
		//std::stringstream ss;
		//ss << path;
		//ss.fill('0');
		//ss.width(5);
		////if(surfNumId%10!=0) return;
		//ss << surfNumId;
		//surfNumId++;
		////if(surfNumId%3!=0) return;
		//ss << ".obj";
		//std::string file_path = ss.str();
		string file_path = path;
		ofstream outSurf(file_path);
		map<int, int> meshToSurf;//(meshTemp.surfVerts.size());
		for (int i = 0; i < meshTemp.surfVerts.size(); i++) {
			const auto& pos = meshTemp.vertices[meshTemp.surfVerts[i]];
			outSurf << "v " << pos.x() << " " << pos.y() << " " << pos.z() << endl;
			meshToSurf[meshTemp.surfVerts[i]] = i;
		}

		for (int i = 0; i < meshTemp.surface.size(); i++) {
			const auto& tri = meshTemp.surface[i];
			outSurf << "f " << meshToSurf[tri[0]] + 1 << " " << meshToSurf[tri[1]] + 1 << " " << meshToSurf[tri[2]] + 1 << endl;
		}
		outSurf.close();
	}

	class mesh_obj {
	public:
		vector<Vector3d> vertexes;
		vector<Vector3d> normals;
		vector<Vector3i> facenormals;
		vector<Vector3i> faces;
		int vertexNum;
		int faceNum;
		void InitMesh(int type, double scale);
		bool load_mesh(const ::std::string& filename, double scale);
	};

	class model_obj {
	public:
		vector<mesh_obj> meshes;
		vector<string> names;
		bool load_model(const ::std::string& filename);
	};

	class model_tet {
	public:
		vector<mesh3D> mesh3Ds;
		//vector<vector<Vector4i>> surfaces;
		vector<string> names;
		bool load_model(const ::std::string& filename, int offset);
		void calculate_surface();
	};

	class fiber_obj {
	public:
		vector<mesh_obj> muscles;
		vector<mesh_obj> tendonIns;
		vector<mesh_obj> tendonOuts;
		//vector<string> names;
		bool load_model(const ::std::string* filename);
	};

}

#endif // !FEM_MESH.H