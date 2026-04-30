#pragma once

#include "bsfwd.h"

#include "BlockAssembleUtils.h"
#include "BSTargetInfo.h"
#include "SolverConfig.h"
#include "Cache.h"

#include "Geometry/Primitives.h"

#include "Energy/Energy.h"

#include "Utility/Cholmod/CholmodSolver.h"

#include "IPC/mesh.h"
#include "IPC/collisionUtil.h"
#include "IPC/IPCdistanceFuncs.h"

// undef conflict macros
#if not defined BSIPC_DISABLE_RENDERER
#undef CreateWindow
#include "Altrar/src/Renderer.h"
#include "Altrar/src/Geometry/Mesh.h"
#endif

#define BSIPC_USE_EIGEN_SOLVER

// check sanity of defines
#if defined BSIPC_USE_EIGEN_SOLVER && defined BSIPC_USE_SUITESPARSE_SOLVER
#error "Both Eigen and SuiteSparse solvers are defined. Please choose one."
#endif

namespace BSIPC
{
    // TODO change all pass-by-value `DMat`s to be members of `Solver` class
    class Solver
    {
    public:
        Solver(SolverConfig config);

        /* ***************************** Contact Barrier ***************************** */

        // These tests highly rely on the interaction between Linear/BS surfaces.
        // For convenience they are not moved to [Energy] folder
#if defined BSIPC_NUMERIC_TEST
        DMat NumericalBarrierGrad_Trig();
        DMat NumericalBarrierHess_Trig();

        DMat NumericalBarrierGrad_BS();
        DMat NumericalBarrierHess_BS();
#endif

        /* ***************************** Solver non-PD Helper ***************************** */

        void FillGlobalGradient(const DMat& localGrad, BSIPC_OUT DMat& globalGrad, const std::vector<UInt> indices, tbb::mutex& mutex) const;
        
        // Returns the assembled hessian multiplied by timestep squared
        void FillElasticityHessTimestepSq(
            const DMat& localHess, BSIPC_OUT SpMatData& globalHessData, const std::vector<UInt> indices, UInt quadPointIndex
        ) const;

        void FillSeamHess(const DMat& localHess, BSIPC_OUT SpMatData& globalHessData, 
            const std::vector<UInt> indicesI, const std::vector<UInt> indicesJ, UInt bsOffsetI, UInt bsOffsetJ,
            UInt hessEntryOffset) const;

        void FillLinearHess(const tbb::concurrent_vector<SpMatEntry>& trigHessData, BSIPC_OUT SpMatData& globalHessData, UInt offset) const;
        void FillAnimatedLinearHess(
            const tbb::concurrent_vector<SpMatEntry>& bsAnimHessData, const tbb::concurrent_vector<SpMatEntry>& animBsHessData,
            const tbb::concurrent_vector<SpMatEntry>& animAnimHessData, BSIPC_OUT SpMatData& globalHessData, UInt offset
        ) const;

        // TODO for efficiently fill in mass matrix
        // Takes matrix entry index (a, b) starting from 0, returns nonzero entries [preceding] it
        // [preceding] is row major, i.e. either row index < b, or row index = b and column index < a
        //UInt PrecedingEntriesInHess(Int a, Int b) const;

        /* ***************************** BC Treatments ***************************** */

        // Fixes control points. according to SolverConfig
        inline void DBCOnGradient(DMat& gradient);

        // returns true if the entry should be fixed in Hess, i.e., the corresponding control point is in the DBC set
        // [indexI] and [indexJ] are the indices in the global Hessian
        Bool ShouldFixHessEntry(UInt indexI, UInt indexJ) const;

        // Due to special treatments on SpMat, DBC on Hessian is done when assembling the SpMatData.
        // See Solver::FillInGlobalHessian and Solver::InertiaHess for details
        // See Energy::ShouldFixHessEntry

        Bool ActivateAnimatedMesh() const;

        /* ***************************** Incremental Potential ***************************** */

        // Evaluated at hypPos. In terms of notation, hypPos = x_tilde, estPos = x_n + dt * v_n
        Float IPVal(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos);
        DMat IPGrad(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos);    // Should return a 3n-by-1 column vector
        SpMat IPHess(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos);

        // Subroutines of IPHessOptimized
        void _ContactBarrierHess(BSIPC_OUT SpMat& hess);
        void ContactBarrierHessReordered(BSIPC_OUT SpMat& hess);
        void ContactBarrierHessReorderedOpt(BSIPC_OUT SpMat& hess);
        void _ContactBarrierHessTriplet(BSIPC_OUT SpMat& hess);

#if defined BSIPC_NUMERIC_TEST
        void BdHessNumericTest(const SpMatData& bdData);                                  // Hessian Assembly Test

        void SsGradBlockNumericTest(const QuadPoint& quadPoint, std::vector<UInt> indices, DMat blockGrad, const BSTargetInfo& info);
        void SsHessBlockNumericTest(const QuadPoint& quadPoint, std::vector<UInt> indices, DMat blockHess, const BSTargetInfo& info);
        void BdGradBlockNumericTest(const QuadPoint& quadPoint, std::vector<UInt> indices, DMat blockGrad, const BSTargetInfo& info);
        void BdHessBlockNumericTest(const QuadPoint& quadPoint, std::vector<UInt> indices, DMat blockHess, const BSTargetInfo& info);
#endif

#if defined BSIPC_NUMERIC_TEST
        void BarrierGradHessNumericTest();
        DMat NumericalIPGrad(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos);

        // When conducting [NumericalIPHess] test, first switch Ss/Bd Hess to non-SPD version.
        DMat NumericalIPHess(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos);
#endif

        /* ***************************** Utils ***************************** */

        void InitSimulation();

        void UpdateContactStiffness();

        // Returns false if the simulation needs to continue, and true if the simulation has reached pre-designated max steps
        Bool StepForward();                                             
        void SolveSubIP(BSIPC_OUT std::vector<Vec3>& curTestBasis, const std::vector<Vec3>& estPos);

        void SolveKKTSystemSchur(const SpMat& hess, const DMat& grad, BSIPC_OUT DMat& x) const;
        Bool SolveKKTSystemWithDofReduction(const SpMat& hess, const DMat& grad, BSIPC_OUT DMat& x) const;
        void PropagateSearchDirToTrigMesh(const DMat& searchDir);
        void PropagateVelToTrigMesh(const std::vector<Vec3>& vel);
        DMat AccumulateTrigGrad(const std::vector<Vec3>& trigGrad) const;

        // Writing to log.txt
        inline void LogSep() const;
        void AppendToLog(String) const;
        void WriteToLog() const;

        // Writing to step_{frame number}.yaml
        // TODO rewrite this part for multiple BSSurfaces
        void WriteToStepLog(const String filename) const;

        // Need to read step log again as in SolverConfig::LoadExperimentLayout the vels/pos of nodes are not processed
        void ReadFromStepLog(const String path);
        void ReadFromStepLog();

        // Called when setting up simulation
        void InitMesh();                                                // Since the indices are fixed, they can be initialized here
        void InitAnimatedMesh();                                        // Process animated meshes. Need to be called after [InitMesh]

        // Requires the mass matrix to be initialized
        void InitContactMesh(std::vector<std::pair<UInt, UInt>> bsTrigSizes);

        // Optimizations for Hessian assembly
        void InitHessCache();
        void PrecalculateBendingHess();
        void PrecomputeKKTComponents();
        void PreAssembleElasticityHessTemplate();

        // Add seam constraints to filtered contact set
        void FilterSeamContacts();

        // Convert seam points from canonical coordiantes to our coordinates
        void NormalizeSeamCoordinates();

        // Called each Newton iteration, when the underlying B-spline surface needs to be updated
        void UpdateControlPoints(const std::vector<Vec3>& newPos);          // Update the mesh with new positions
        void UpdateAnimatedMesh(const std::vector<Vec3>& pos);
        void UpdateMovingDBC(std::vector<Vec3>& newPos, UInt curStep);
        void UpdateTrigMesh(BSIPC_OUT IPC::mesh3D& mesh) const;
        void UpdateContactMesh();
        void UpdatePatchOrderMesh();
        void UpdateContactActiveSets();

        // Construct triangle indices in patch order. Should only be used when translating yaml to obj
        void InitPatchOrderMesh(const std::vector<std::pair<UInt, UInt>>& bsTrigSizes);

#if not defined BSIPC_DISABLE_RENDERER
        void InitRenderMesh(std::vector<std::pair<UInt, UInt>> bsTrigSizes);
        
        // Called each step, for demo rendering
        void UpdateRenderMesh();                                        // So when updating mesh, only the vertex positions need to be updated

        inline const ATR::Mesh& GetMesh() const { return renderMesh; }
#endif

        /* ***************************** Getter/Setter Utils ***************************** */
#if not defined BSIPC_DISABLE_RENDERER
        inline void BindRenderer(ATR::Renderer* renderer) { this->renderer = renderer; }
#endif

        inline const SolverConfig& GetConfig() const { return config; }

        inline const TriangularMesh& GetPatchOrderMesh(UInt idx) const { return patchOrderMesh[idx]; }
        inline const std::vector<Vec3>& GetControlPoints(UInt index) const { return this->GetTarget(index)->GetControlPoints(); }

        // Wrapper for getting underlying BS surfaces
        inline UInt const GetTargetCnt() { return static_cast<UInt>(this->targets.size()); }
        inline BSSurface* const GetTarget(UInt index) { return this->targets[index].target; }
        inline const BSSurface* const GetTarget(UInt index) const { return this->targets[index].target; }
        inline const std::vector<BSTargetInfo>& GetTargetsInfo() const { return this->targets; }

        inline void AddTarget(BSSurface* target, Vec3 initVel)
        { 
            this->targets.emplace_back(target, static_cast<UInt>(this->targets.size()));
            UInt nodeCnt = target->GetControlPointCnt();
            this->bsVertCntSum += nodeCnt;
            this->bsPatchCntSum += target->GetPatchCnt();

            for (UInt i = 0; i != nodeCnt; ++i)
                this->stepCache.estVels.push_back(initVel);
        }

        inline void AddTarget(BSSurface* target, const std::vector<Vec3>& initVel)
        {
            this->targets.emplace_back(target, static_cast<UInt>(this->targets.size()));
            UInt nodeCnt = target->GetControlPointCnt();
            this->bsVertCntSum += nodeCnt;
            this->bsPatchCntSum += target->GetPatchCnt();

            BSIPC_ASSERT(initVel.size() == nodeCnt, "specified init velocity has size inconsistent with nodeCnt");

            for (UInt i = 0; i != nodeCnt; ++i)
                this->stepCache.estVels.push_back(initVel[i]);
        }

        inline void FormTrigInIPCMesh(Int v1, Int v2, Int v3);

        // [indI] and [indJ] are (1-dimensional) vertex indices in the contact mesh
        inline UInt GetContactMeshIndex(UInt indI, UInt indJ, const BSTargetInfo& targetInfo) const;

        inline UInt CurStep() const { return this->stepCache.CurStep(); }

    private:
        SolverConfig config;
        std::vector<BSTargetInfo> targets;
        std::vector<FixedMeshInfo> fixedLinearMeshes;
        StepCache stepCache;

        // TODO may want to optimize by precomputing quantities at quadrature points
        QuadratureCache quadCache;
        std::vector<Float> massMatDiagEntries;

        std::vector<TriangularMesh> patchOrderMesh;

        std::vector<Mat<27, 1>> localElasticityGrad;                    // Local 27 * 1 elasticity gradient
        std::vector<Mat<27, 27>> localBendingHess;                      // Precalculate constant bending hessian for each patch
        std::vector<Mat<27, 27>> localElasticityHess;                   // Local 27 * 27 elasticity Hessian

        // KKT system components [See PrecomputeKKTComponents()]
        UInt seamConstraintCnt;
        SpMat seamConstraintCoeffs;                                     // The coefficient matrix C s.t. C.x = 0
        SpMat seamConstraintCoeffsExtended;                             // Similar to above, but expanded to the size of Hessian including animation sequence DoF

        // KKT system projections [See PrecomputeKKTComponents()]
        SpMat rowReorderingMat;                                         // The R matrix in derivation
        SpMat rowReorderingMatExtended;                                 // Similar to above, but expanded to the size of Hessian including animation sequence DoF
        SpMat constrainedDofProjMat;                                    // The Z matrix in derivation. Z is restricted to only constrained DoFs, and is therefore dense
        SpMat constrainedDofProjMatExtended;                            // Similar to above, but expanded to the size of Hessian including animation sequence DoF

        //std::vector<std::vector<ElasticityGradBlockSource>> elasticityGradSources; // map: gradient entry index -> list of block gradient sources
        std::vector<std::vector<ElasticityHessEntrySource>> elasticityHessSources; // map: hessian [data] array index -> list of block hessian sources
        SpMat elasticityHess;                                           // Assembled global elasticity Hessian. Modified based on template.
        SpMat elasticityHessTemplate;                                   // Precomputed, w/ full sparsity pattern and mass matrix value; should not be modified during runtime.

        // IPC/Collision interfaces

        /// The sources recorded in the following only accounts for B-spline surfaces. DoF of animation sequence is not included.
        std::vector<std::vector<SourceInfo>> linearVertexSources;       // map: linear mesh vertex index -> list of BS controls point index the vertex is supported by; and its weight
        std::vector<std::vector<SourceInfo>> bsVertexSupports;          // map: BS control point index -> list of linear mesh vertices in its support; and its weight

        UInt bsVertCntSum;
        UInt trigVertCntSum;
        UInt bsPatchCntSum;

        Float contactStiffness;

        UInt inertiaHessEntryCnt, potHessEntryCnt;

        UInt curIndexInHessCache;                                       // records which entry should be recorded next iteration
        std::vector<SpMatEntry> hessEntries;

        // NOTE: IPC library integrates all contact tests in one mesh, so ths number of vertices will be larger than those belonging to BS surfaces
        //  so does the velocity/searchDir vectors
        IPC::mesh3D contactMesh;
        IPC::Ground ground;
        IPC::SpatialHash spatialHash;

        std::vector<Vec3> contactMeshBarrierGrad;
        IPC::BHessian contactMeshBarrierHess;

        std::vector<Vec3> contactMeshSearchDir;

        // mesh animation sequence
        std::optional<AnimatedMeshCache> animatedMeshCache;
        std::vector<IPC::mesh3D> animatedMeshPoses;
        Bool includeAnimatedMeshInSystem = false;
            
        // Utils
        UInt simlTime;

#if not defined BSIPC_DISABLE_RENDERER
        // Rendering
        ATR::Renderer* renderer;
        ATR::Mesh renderMesh;                             // Mesh for rendering. For simplicity of optimization this is separate from the triangular mesh for contact
#endif
    };

#if not defined BSIPC_DISABLE_RENDERER
    inline ATR::Vec3 TransformVec3(const Vec3& vec) { return ATR::Vec3(vec[0], vec[1], vec[2]); }
#endif

    inline void Solver::LogSep() const
    {
        this->AppendToLog("========================================\n");
    }

    inline void Solver::FormTrigInIPCMesh(Int v1, Int v2, Int v3)
    {
        this->contactMesh.triangles.emplace_back(v1, v2, v3);
    }

    inline UInt Solver::GetContactMeshIndex(UInt indI, UInt indJ, const BSTargetInfo& info) const
    {
        BSIPC_ASSERT(info.ValidateTrig(), "Contact mesh offset not set.");

        UInt bsWidth = info.target->GetResolutionX();
        UInt perPatch = this->config.perPatchMesh;
        UInt resolutionX = (bsWidth - 1) * perPatch + 1;
        UInt index = indI + indJ * resolutionX;

        index += info.contactTrigVertIndexOffset.value();

        return index;
    }
}
