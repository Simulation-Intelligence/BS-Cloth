#pragma once

#include "bsfwd.h"

#include "fkyaml/node.hpp"
#include "Geometry/Transformation.h"

namespace BSIPC
{
    enum class QuadScheme
    {
        LOCAL_UNIFORM,
        LOCAL_EDGE_DENSE,
        LOCAL_CHESSBOARD,
        LOCAL_CHESSBOARD_DENSE,

        LOCAL_DUAL_CHESSBOARD,

        GLOBAL_UNIFORM,
        GLOBAL_LIMITED,
        GLOBAL_EDGE_DENSE,
        GLOBAL_CHESSBOARD
    };

    enum class FixedMeshType
    {
        TETRAHEDRAL,
        TRIANGLE
    };

    QuadScheme QuadStrToScheme(const String& quadStr);
    String QuadSchemeToStr(QuadScheme scheme);

    struct MovingNode
    {
    public:
        UInt index;
        Vec3 velocity;
        UInt untilFrame;

        Bool operator==(const MovingNode& other) const { return index == other.index; }
        Bool operator==(UInt index) const { return index == index; }
    };

    enum StartConfig
    {
        VERTICAL,
        FLAT
    };

    struct FixedMeshInfo
    {
        FixedMeshType type;
        String meshName;
        std::vector<LinearTransformation> transformations;
    };

    using SeamPoint = std::pair<Vec2, Vec2>;

    struct SeamInfo
    {
        std::array<UInt, 2> seamSurfIndices = { 0, 0 };
        std::vector<SeamPoint> seamPoints;
        String seamingPointsFileName;
        Bool inCanonicalCoord = false;
    };

    struct AnimatedMeshInfo
    {
        String dirName;
        String pattern;
        UInt length;
        Float stiffness;
        Float tolerance;
        std::vector<LinearTransformation> transformations;
    };

    struct BSSurfaceInfo
    {
        UInt idx;

        UInt resolutionX;
        UInt resolutionY;

        UInt linearResolutionX;
        UInt linearResolutionY;

        std::vector<LinearTransformation> transformations;

        Float restSizeX;
        Float restSizeY;

        Float initRatioX;
        Float initRatioY;

        Vec3 initVel;
        std::vector<Vec3> initVelocities;
        Vec3 restPos;

        std::unordered_set<UInt> fixedNodes;
        std::vector<MovingNode> movingNodes;

        Bool irregularDomain;
        String restNodeFileName;
        std::vector<Vec2> restNodePos;
        std::vector<Vec3> startNodePos;

        StartConfig startConfig;

        Float volume;
        Float surfaceArea;

        BSSurfaceInfo() :
            idx(std::numeric_limits<UInt>::max()),
            resolutionX(20), resolutionY(20), linearResolutionX(20), linearResolutionY(40),
            restSizeX(0.8), restSizeY(0.8),
            initRatioX(1.0), initRatioY(1.0), initVel(Vec3::Zero()), irregularDomain(false),
            startConfig(StartConfig::FLAT), volume(0.00064), surfaceArea(0.64) {}
    };

    struct SolverConfig
    {
        Bool continueOnStep;
        UInt startStep;
        fkyaml::node configRootNode;

        Float gravConst;

        Float ssYoungsModulus;
        Float bdYoungsModulus;
        Float poissonRatio;

        Float bu, bv;               // The shrink/stretch parameters
        Float stretchStiffness;
        Float shearStiffness;
        Float strainRate;
        Float bendingStiffness;
        Float mu, lambda;
        std::optional<Float> contactStiffness;
        Float frictionMu;
        Float seamStrength;

        Bool enableContact;
        Bool disableInterpatchContact;
        Bool writeIterMesh;
        Float contactThreshold;

        Float groundHeight;

        Float timestep;
        Float tolerance;

        UInt maxLSIterations;
        UInt maxNewtonIterations;

        QuadScheme quadScheme;
        UInt quadOrder;

        UInt maxSteps;

        Float thickness;

        Float density;

        std::vector<BSSurfaceInfo> surfacesInfo;
        std::vector<FixedMeshInfo> fixedMeshInfos;
        std::vector<SeamInfo> seamInfos;
        std::optional<AnimatedMeshInfo> animatedMeshInfo;

        UInt perPatchMesh;

        Bool screenshot;
        Bool writeIPCMesh;
        Bool cacheQuadPoints;
        Bool useCachedQuadPoints;

        static SolverConfig DefaultInfo();

        void LoadRestNodePos(const String& path, BSIPC_OUT BSSurfaceInfo& info) const;
        void LoadStartNodePos(const String& path, BSIPC_OUT BSSurfaceInfo& info) const;
        void LoadInitVelocities(const String& path, BSIPC_OUT BSSurfaceInfo& info) const;
        void LoadSurfaceInfos(const fkyaml::node& surfacesInfoNode, BSIPC_OUT std::vector<BSSurfaceInfo>& surfacesInfo) const;
        void LoadSeamInfo(const String& path, BSIPC_OUT SeamInfo& info) const;
        void LoadFixedMeshInfo(const fkyaml::node& parent, BSIPC_OUT std::vector<FixedMeshInfo>& info) const;
        void LoadAnimatedMeshInfo(const fkyaml::node& parent, BSIPC_OUT std::optional<AnimatedMeshInfo>& info) const;

        void LoadExperimentLayout(const String& path);
        void LoadExperimentLayout();
        void SerializeConfig();

        void LoadTransformations(const fkyaml::node& parent, BSIPC_OUT std::vector<LinearTransformation>& transformations) const;
        void SerializeTransformations(BSIPC_OUT fkyaml::node& parent, const std::vector<LinearTransformation>& transformations) const;

        std::optional<Vec3> MovingNodeVelocity(UInt index, UInt nodeIndex) const;
    };

    namespace Config
    {
        // Applicable for all configurations
        inline const char* CONTINUE_ON_STEP_STR = "continue_on_step";
        inline const char* STEP_PATH_STR = "step_path";
        inline const char* DENSITY_STR = "density";
        inline const char* MAX_STEPS_STR = "max_steps";
        inline const char* MAX_LS_ITERATIONS_STR = "max_ls_iterations";
        inline const char* MAX_NEWTON_ITERATIONS_STR = "max_newton_iterations";

        inline const char* BS_SURFACES_STR = "bs_surfaces";
        inline const char* IDX_STR = "idx";
        inline const char* RESOLUTION_STR = "resolution";
        inline const char* LINEAR_RESOLUTION_STR = "linear_resolution";
        inline const char* X_STR = "x";
        inline const char* Y_STR = "y";
        inline const char* REST_SIZE_STR = "rest_size";
        inline const char* INIT_RATIO_STR = "init_ratio";
        inline const char* START_CONFIG_STR = "start_config";
        inline const char* FLAT_STR = "flat";
        inline const char* VERTICAL_STR = "vertical";
        inline const char* FIXED_NODES_STR = "fixed_nodes";
        inline const char* MOVING_NODES_STR = "moving_nodes";
        inline const char* MOVING_NODE_ID_STR = "id";
        inline const char* MOVING_NODE_VELOCITY_STR = "vel";
        inline const char* MOVING_NODE_UNTIL_FRAME = "until_frame";
        inline const char* REST_NODE_FILE_STR = "rest_node_file";
        inline const char* START_NODE_FILE_STR = "start_node_file";
        inline const char* INIT_VELOCITY_STR = "init_velocity";
        inline const char* INIT_VELOCITY_FILE_STR = "init_velocity_file";
        inline const char* REST_POS_STR = "rest_pos";

        inline const char* SS_QUAD_POINTS_STR = "ss_quad_points";
        inline const char* BD_QUAD_POINTS_STR = "bd_quad_points";
        inline const char* MASS_MATRIX_STR = "mass_matrix";
        inline const char* PATCH_START_INDEX_STR = "patch_start_index";
        inline const char* USE_CACHED_QUAD_POINTS_STR = "use_cached_quad_points";

        inline const char* FIXED_TET_MESHES_STR = "fixed_tet_meshes";
        inline const char* FIXED_TRIG_MESHES_STR = "fixed_trig_meshes";
        inline const char* NAME_STR = "name";
        inline const char* TRANSFORMATIONS_STR = "transformations";
        inline const char* TRANSFORMATION_TYPE_STR = "type";
        inline const char* TRANSFORMATION_CONTENT_STR = "content";
        inline const char* TRANSLATION_STR = "translation";
        inline const char* ROTATION_STR = "rotation";
        inline const char* SCALING_STR = "scaling";

        inline const char* ANIMATED_MESHES_STR = "animated_meshes";
        inline const char* DIRECTORY_STR = "directory";
        inline const char* PATTERN_STR = "pattern";
        inline const char* LENGTH_STR = "length";
        inline const char* STIFFNESS_STR = "stiffness";
        inline const char* DIFF_STR = "diff";
        inline const char* VERTICES_STR = "vertices";

        inline const char* QUAD_SCHEME_STR = "quad_scheme";
        inline const char* QUAD_LOCAL_STR = "local";
        inline const char* QUAD_LOCAL_EDGE_DENSE_STR = "local-edge-dense";              // edge 3-3, interior 2-2
        inline const char* QUAD_LOCAL_CHESSBOARD_STR = "local-chessboard";              // edge 3-3, interior 1-2/2-1
        inline const char* QUAD_LOCAL_CHESSBOARD_DENSE_STR = "local-chessboard-dense";  // edge 3-3, interior 2-2/1-1
        inline const char* QUAD_LOCAL_DUAL_CHESSBOARD_STR = "local-dual-chessboard";    // edge 1-1, interior 1-2/2-1 int the dual grid
        inline const char* QUAD_GLOBAL_STR = "global";
        inline const char* QUAD_GLOBAL_LIMITED_STR = "global-limited";
        inline const char* QUAD_GLOBAL_EDGE_DENSE_STR = "global-edge-dense";            // edge 3-3, interior 2-2
        inline const char* QUAD_GLOBAL_CHESSBOARD_STR = "global-chessboard";            // edge 3-3, interior 1-2/2-1
        inline const char* QUAD_ORDER_STR = "quad_order";

        inline const char* SEAMS_STR = "seams";
        inline const char* IDX_PAIR_STR = "idx_pair";
        inline const char* SEAMING_POINTS_FILE_STR = "seaming_points_file";
        inline const char* IN_CANONICAL_COORD_STR = "in_canonical_coord";

        inline const char* PER_PATCH_MESH_STR = "per_patch_mesh";
        inline const char* SCREENSHOT_STR = "screenshot";
        inline const char* WRITE_IPC_MESH_STR = "write_ipc_mesh";
        inline const char* WRITE_ITER_MESH_STR = "write_iter_mesh";
        inline const char* CACHE_QUAD_POINTS_STR = "cache_quad_points";
        inline const char* ENABLE_CONTACT_STR = "enable_contact";
        inline const char* CONTACT_THRESHOLD_STR = "contact_threshold";
        inline const char* GROUND_HEIGHT_STR = "ground_height";
        inline const char* THICKNESS_STR = "thickness";
        inline const char* GRAVITY_STR = "gravity";
        inline const char* DISABLE_INTERPATCH_CONTACT_STR = "disable_interpatch_contact";

        inline const char* STRETCHING_YOUNGS_MODULUS_STR = "stretching_youngs_modulus";
        inline const char* BENDING_YOUNGS_MODULUS_STR = "bending_youngs_modulus";
        inline const char* SHEAR_STRECTH_RATIO_STR = "shear_stretch_ratio";
        inline const char* POISSON_RATIO_STR = "poisson_ratio";
        inline const char* CONTACT_STIFFNESS_STR = "contact_stiffness";
        inline const char* FRICTION_MU_STR = "friction_mu";
        inline const char* SEAM_STRENGTH_STR = "seam_strength";

        inline const char* TIMESTEP_STR = "timestep";
        inline const char* RESIDUE_TOLERANCE_STR = "residue_tolerance";

        // Read only if continue_on_step is true
        inline const char* CONTROL_NODES_STR = "control_nodes";
        inline const char* CURRENT_STEP_STR = "current_step";
        inline const char* EST_VELOCITY_STR = "est_velocity";
    };

    inline bool IsQuadSchemeGlobal(QuadScheme scheme)
    {
        return scheme == QuadScheme::GLOBAL_UNIFORM ||
            scheme == QuadScheme::GLOBAL_LIMITED ||
            scheme == QuadScheme::GLOBAL_EDGE_DENSE ||
            scheme == QuadScheme::GLOBAL_CHESSBOARD;
    }

    inline bool IsQuadSchemeLocal(QuadScheme scheme)
    {
        return scheme == QuadScheme::LOCAL_UNIFORM ||
            scheme == QuadScheme::LOCAL_EDGE_DENSE ||
            scheme == QuadScheme::LOCAL_CHESSBOARD ||
            scheme == QuadScheme::LOCAL_CHESSBOARD_DENSE;
    }
}
