#include "bspch.h"

#include "SolverConfig.h"
#include "Utility/Loader/YAMLLoader.h"

namespace BSIPC
{
    QuadScheme QuadStrToScheme(const String& quadStr)
    {
        if (quadStr == Config::QUAD_LOCAL_STR)
            return QuadScheme::LOCAL_UNIFORM;
        else if (quadStr == Config::QUAD_LOCAL_EDGE_DENSE_STR)
            return QuadScheme::LOCAL_EDGE_DENSE;
        else if (quadStr == Config::QUAD_LOCAL_CHESSBOARD_STR)
            return QuadScheme::LOCAL_CHESSBOARD;
        else if (quadStr == Config::QUAD_LOCAL_CHESSBOARD_DENSE_STR)
            return QuadScheme::LOCAL_CHESSBOARD_DENSE;
        else if (quadStr == Config::QUAD_GLOBAL_STR)
            return QuadScheme::GLOBAL_UNIFORM;
        else if (quadStr == Config::QUAD_GLOBAL_LIMITED_STR)
            return QuadScheme::GLOBAL_LIMITED;
        else if (quadStr == Config::QUAD_GLOBAL_EDGE_DENSE_STR)
            return QuadScheme::GLOBAL_EDGE_DENSE;
        else if (quadStr == Config::QUAD_GLOBAL_CHESSBOARD_STR)
            return QuadScheme::GLOBAL_CHESSBOARD;
        else if (quadStr == Config::QUAD_LOCAL_DUAL_CHESSBOARD_STR)
            return QuadScheme::LOCAL_DUAL_CHESSBOARD;
        else
        {
            BSIPC_WARN("Invalid quad_scheme in config, using default: global uniform");
            return QuadScheme::GLOBAL_UNIFORM;
        }
    }

    String QuadSchemeToStr(QuadScheme scheme)
    {
        if (scheme == QuadScheme::LOCAL_UNIFORM)
            return Config::QUAD_LOCAL_STR;
        else if (scheme == QuadScheme::LOCAL_EDGE_DENSE)
            return Config::QUAD_LOCAL_EDGE_DENSE_STR;
        else if (scheme == QuadScheme::LOCAL_CHESSBOARD)
            return Config::QUAD_LOCAL_CHESSBOARD_STR;
        else if (scheme == QuadScheme::LOCAL_CHESSBOARD_DENSE)
            return Config::QUAD_LOCAL_CHESSBOARD_DENSE_STR;
        else if (scheme == QuadScheme::LOCAL_DUAL_CHESSBOARD)
            return Config::QUAD_LOCAL_DUAL_CHESSBOARD_STR;
        else if (scheme == QuadScheme::GLOBAL_UNIFORM)
            return Config::QUAD_GLOBAL_STR;
        else if (scheme == QuadScheme::GLOBAL_LIMITED)
            return Config::QUAD_GLOBAL_LIMITED_STR;
        else if (scheme == QuadScheme::GLOBAL_EDGE_DENSE)
            return Config::QUAD_GLOBAL_EDGE_DENSE_STR;
        else if (scheme == QuadScheme::GLOBAL_CHESSBOARD)
            return Config::QUAD_GLOBAL_CHESSBOARD_STR;
        else
        {
            BSIPC_WARN("Invalid quad_scheme in config, using default: global uniform");
            return Config::QUAD_GLOBAL_STR;
        }
    }

    SolverConfig SolverConfig::DefaultInfo()
    {
        Float ssYoungsModulus = 1e5;
        Float poissonRatio = 0.45;

        Float stretchStiffness = ssYoungsModulus / (2 * (1 + poissonRatio));
        Float shearStiffness = 0.3 * stretchStiffness;
        const Float thickness = 1e-3;
        Float bdYoungsModulus = 1e6;
        Float bendingStiffness = bdYoungsModulus * thickness * thickness * thickness / (24 * (1 - poissonRatio * poissonRatio));

        return SolverConfig{
            .continueOnStep = false,
            .startStep = 0,
            .configRootNode = fkyaml::node(fkyaml::node_type::MAPPING),
            .gravConst = 9.81,
            .ssYoungsModulus = ssYoungsModulus,
            .bdYoungsModulus = bdYoungsModulus,
            .poissonRatio = poissonRatio,
            .bu = 1.,
            .bv = 1.,
            .stretchStiffness = stretchStiffness,
            .shearStiffness = shearStiffness,
            .strainRate = 0.,
            .bendingStiffness = bendingStiffness,
            .frictionMu = 0.1,
            .seamStrength = 1.,
            .enableContact = true,
            .disableInterpatchContact = false,
            .contactThreshold = 1e-3,
            .groundHeight = -5.0,
            .timestep = 1e-2,
            .tolerance = 1e-2,
            .maxLSIterations = 500,
            .quadScheme = QuadScheme::LOCAL_EDGE_DENSE,
            .quadOrder = 6,
            .maxSteps = 1000,
            .thickness = 1e-3,
            .density = 200.,
            .surfacesInfo = {},
            .perPatchMesh = 2,
            .screenshot = false,
            .writeIPCMesh = false,
            .cacheQuadPoints = false,
            .useCachedQuadPoints = false,
        };
    }

    void SolverConfig::LoadRestNodePos(const String& path, BSIPC_OUT BSSurfaceInfo& info) const
    {
        std::ifstream ifs(path);

        try
        {
            if (!ifs)
            {
                std::string msg = "error opening config file " + path;
                throw fkyaml::exception(msg.c_str());
            }
            fkyaml::node root = fkyaml::node::deserialize(ifs);

            BSIPC_ASSERT(root.is_sequence(), "node rest position file should be a sequence of 2-element arrays");

            std::vector<Vec2> restNodePos;
            for (auto& node : root)
            {
                BSIPC_ASSERT(node.is_sequence() && node.size() == 2, 
                    "node rest position should be a sequence of 2-element arrays");
                
                Vec2 restPos = Vec2(node[0].get_value<Float>(), node[1].get_value<Float>());

                restNodePos.push_back(restPos);
            }

            BSIPC_ASSERT(
                restNodePos.size() == info.resolutionX * info.resolutionY,
                fmt::format("node rest position size ({}) does not match the resolution ({} x {})", restNodePos.size(), info.resolutionX, info.resolutionY
            ));

            info.restNodePos = restNodePos;
        }
        catch (fkyaml::exception& e)
        {
            BSIPC_WARN(e.what());
            BSIPC_INFO("Error loading rest nodes position from file \"{}\". Ignoring.", path);
        }
    }

    void SolverConfig::LoadStartNodePos(const String& path, BSSurfaceInfo& info) const
    {
        std::ifstream ifs(path);

        try
        {
            if (!ifs)
            {
                std::string msg = "error opening config file " + path;
                throw fkyaml::exception(msg.c_str());
            }
            fkyaml::node root = fkyaml::node::deserialize(ifs);

            BSIPC_ASSERT(root.is_sequence(), "node rest position file should be a sequence of 3-element arrays");

            std::vector<Vec3> startNodePos;
            for (auto& node : root)
            {
                BSIPC_ASSERT(node.is_sequence() && node.size() == 3,
                    "node rest position should be a sequence of 3-element arrays");

                Vec3 startPos = Vec3(node[0].get_value<Float>(), node[1].get_value<Float>(), node[2].get_value<Float>());

                startNodePos.push_back(startPos);
            }

            BSIPC_ASSERT(startNodePos.size() == info.resolutionX * info.resolutionY,
                "node rest position size does not match the resolution");

            info.startNodePos = startNodePos;
        }
        catch (fkyaml::exception& e)
        {
            BSIPC_WARN(e.what());
            BSIPC_INFO("Error loading rest nodes position from file \"{}\". Ignoring.", path);
        }
    }

    void SolverConfig::LoadInitVelocities(const String& path, BSIPC_OUT BSSurfaceInfo& info) const
    {
        std::ifstream ifs(path);

        try
        {
            if (!ifs)
            {
                std::string msg = "error opening config file " + path;
                throw fkyaml::exception(msg.c_str());
            }
            fkyaml::node root = fkyaml::node::deserialize(ifs);

            BSIPC_ASSERT(root.is_sequence(), "node rest position file should be a sequence of 3-element arrays");

            std::vector<Vec3> initVelocities;
            for (auto& node : root)
            {
                BSIPC_ASSERT(node.is_sequence() && node.size() == 3,
                    "node rest position should be a sequence of 3-element arrays");

                Vec3 curInitVelocity = Vec3(node[0].get_value<Float>(), node[1].get_value<Float>(), node[2].get_value<Float>());

                initVelocities.push_back(curInitVelocity);
            }

            info.initVelocities = initVelocities;
        }
        catch (fkyaml::exception& e)
        {
            BSIPC_WARN(e.what());
            BSIPC_INFO("Error loading rest nodes position from file \"{}\". Ignoring.", path);
        }
    }

    void SolverConfig::LoadSurfaceInfos(const fkyaml::node& surfacesInfoNode, BSIPC_OUT std::vector<BSSurfaceInfo>& surfacesInfo) const
    {
        for (auto& surface : surfacesInfoNode)
        {
            BSIPC_ASSERT(surface.is_mapping(), "subnodes of \"bs_surfaces\" in config is not mapping");
            BSSurfaceInfo info;

            //LoadDataFromYaml(info.idx, surface, Config::IDX_STR);

            auto resolutionNode = surface;
            LoadNodeFromYaml(resolutionNode, surface, Config::RESOLUTION_STR);

            LoadDataFromYamlOptional(info.resolutionX, resolutionNode, Config::X_STR);
            LoadDataFromYamlOptional(info.resolutionY, resolutionNode, Config::Y_STR);

            auto linearResolutionNode = surface;
            LoadNodeFromYaml(linearResolutionNode, surface, Config::LINEAR_RESOLUTION_STR);

            LoadDataFromYamlOptional(info.linearResolutionX, linearResolutionNode, Config::X_STR);
            LoadDataFromYamlOptional(info.linearResolutionY, linearResolutionNode, Config::Y_STR);

            LoadDataFromYamlOptional(info.idx, surface, Config::IDX_STR);

            this->LoadTransformations(surface, info.transformations);

            String initVelFileStr;
            Bool specifiedInitVelFile = LoadDataFromYamlOptional(initVelFileStr, surface, Config::INIT_VELOCITY_FILE_STR);
            if (specifiedInitVelFile)
            {
                String inDirInitVelFile = fmt::format("models/{}.yaml", initVelFileStr);
                this->LoadInitVelocities(inDirInitVelFile, info);
            }
            else
            {
                auto initVelNode = surface;
                Bool specifiedInitVel = LoadNodeFromYamlOptional(initVelNode, surface, Config::INIT_VELOCITY_STR);
                if (specifiedInitVel)
                {
                    BSIPC_ASSERT(initVelNode.is_sequence() && initVelNode.size() == 3,
                        "\"init_velocity\" node in config is not a sequence of 3 elements");
                    Vec3 initVel;
                    for (UInt i = 0; i < 3; ++i)
                    {
                        BSIPC_ASSERT(initVelNode[i].is_float_number(), "\"init_velocity\" node has non-floating number data");
                        initVel[i] = initVelNode[i].get_value<Float>();
                    }
                    info.initVel = initVel;
                }
                else
                    info.initVel = Vec3::Zero();
            }

            auto restPosNode = surface;
            Bool specifiedRestPos = LoadNodeFromYamlOptional(restPosNode, surface, Config::REST_POS_STR);
            if (specifiedRestPos)
            {
                BSIPC_ASSERT(restPosNode.is_sequence() && restPosNode.size() == 3,
                    "\"rest_pos\" node in config is not a sequence of 3 elements");
                Vec3 restPos;
                for (UInt i = 0; i < 3; ++i)
                {
                    BSIPC_ASSERT(restPosNode[i].is_float_number(), "\"rest_pos\" node has non-floating number data");
                    restPos[i] = restPosNode[i].get_value<Float>();
                }
                info.restPos = restPos;
            }
            else
                info.restPos = Vec3::Zero();

            auto restSizeNode = surface;
            Bool specifiedRestSize = LoadNodeFromYamlOptional(restSizeNode, surface, Config::REST_SIZE_STR);
            if (specifiedRestSize)
            {
                LoadDataFromYamlOptional(info.restSizeX, restSizeNode, Config::X_STR);
                LoadDataFromYamlOptional(info.restSizeY, restSizeNode, Config::Y_STR);
            }

            String configStr;
            Bool specifiedStartConfig = LoadDataFromYamlOptional(configStr, surface, Config::START_CONFIG_STR);
            if (specifiedStartConfig)
            {
                if (configStr == Config::FLAT_STR)
                    info.startConfig = StartConfig::FLAT;
                else if (configStr == Config::VERTICAL_STR)
                    info.startConfig = StartConfig::VERTICAL;
                else
                    BSIPC_INFO("Start config \"{}\" is not valid. Use [flat] or [vertical]", configStr);
            }

            auto initRatioNode = surface;
            Bool specifiedInitRatio = LoadNodeFromYamlOptional(initRatioNode, surface, Config::INIT_RATIO_STR);

            if (specifiedInitRatio)
            {
                LoadDataFromYamlOptional(info.initRatioX, initRatioNode, Config::X_STR);
                LoadDataFromYamlOptional(info.initRatioY, initRatioNode, Config::Y_STR);
            }

            // Fixed DBC
            auto fixedSpecNode = surface;
            if (LoadNodeFromYamlOptional(fixedSpecNode, surface, Config::FIXED_NODES_STR))
            {
                info.fixedNodes.clear();
                for (auto& node : fixedSpecNode)
                {
                    BSIPC_ASSERT(node.is_integer(), "\"fixed_nodes\" in config is not integer")
                        info.fixedNodes.insert(node.get_value<UInt>());
                }
            }

            // Moving DBC
            auto movingSpecNode = surface;
            if (LoadNodeFromYamlOptional(movingSpecNode, surface, Config::MOVING_NODES_STR))
            {
                BSIPC_INFO("reading moving nodes config");

                for (auto& node : movingSpecNode)
                {
                    BSIPC_ASSERT(node.is_mapping(), "subnodes of \"moving_nodes\" in config is not mapping");

                    MovingNode curMovingNode;
                    LoadDataFromYaml(curMovingNode.index, node, Config::MOVING_NODE_ID_STR);

                    auto velocityNode = node;
                    LoadNodeFromYaml(velocityNode, node, Config::MOVING_NODE_VELOCITY_STR);
                    BSIPC_ASSERT(velocityNode.is_sequence() && velocityNode.size() == 3, "invalid \"velocity\" node data");
                    Vec3 curVelocity;
                    for (UInt i = 0; i < 3; ++i)
                    {
                        BSIPC_ASSERT(velocityNode[i].is_float_number(), "\"velocity\" node has non-floating number data");
                        curVelocity[i] = velocityNode[i].get_value<Float>();
                    }

                    LoadDataFromYaml(curMovingNode.untilFrame, node, Config::MOVING_NODE_UNTIL_FRAME);

                    curMovingNode.velocity = curVelocity;
                    info.movingNodes.push_back(curMovingNode);
                }
            }

            // Irregular domain
            info.irregularDomain = LoadDataFromYamlOptional(info.restNodeFileName, surface, Config::REST_NODE_FILE_STR);
            BSIPC_PEEK(info.irregularDomain);

            if (info.irregularDomain)
            {
                String inDirRestNodeFile = fmt::format("models/{}.yaml", info.restNodeFileName);
                this->LoadRestNodePos(inDirRestNodeFile, info);

                Float totalArea = 0.;
                for (UInt i = 0; i < info.resolutionX - 1; ++i)
                    for (UInt j = 0; j < info.resolutionY - 1; ++j)
                    {
                        Vec2 v1 = info.restNodePos[j * info.resolutionX + i];
                        Vec2 v2 = info.restNodePos[(j + 1) * info.resolutionX + i];
                        Vec2 v3 = info.restNodePos[j * info.resolutionX + (i + 1)];
                        Vec2 v4 = info.restNodePos[(j + 1) * info.resolutionX + (i + 1)];

                        totalArea += trigArea(v1, v2, v3);
                        totalArea += trigArea(v2, v3, v4);
                    }

                info.surfaceArea = totalArea;
                info.volume = info.surfaceArea * this->thickness;
            }
            else
            {
                info.volume = info.restSizeX * info.restSizeY * this->thickness;
                info.surfaceArea = info.restSizeX * info.restSizeY;
            }

            // Start configuration, not necessarily irregular
            String startNodeFile;
            Bool specifiedStartNode = LoadDataFromYamlOptional(startNodeFile, surface, Config::START_NODE_FILE_STR);
            if (specifiedStartNode)
            {
                String inDirStartNodeFile = fmt::format("models/{}.yaml", startNodeFile);
                this->LoadStartNodePos(inDirStartNodeFile, info);
            }

            surfacesInfo.emplace_back(info);
        }
    }

    void SolverConfig::LoadSeamInfo(const String& path, SeamInfo& info) const
    {
        std::ifstream ifs(path);

        try
        {
            if (!ifs)
            {
                std::string msg = "error opening config file " + path;
                throw fkyaml::exception(msg.c_str());
            }
            fkyaml::node root = fkyaml::node::deserialize(ifs);


            LoadDataFromYaml(info.inCanonicalCoord, root, Config::IN_CANONICAL_COORD_STR);

            fkyaml::node seamNode = root;
            LoadNodeFromYaml(seamNode, root, Config::SEAMS_STR);
            BSIPC_ASSERT(seamNode.is_sequence(), "seaming points file should be a sequence of 2-element arrays");

            std::vector<SeamPoint> seamNodes;

            for (auto& node : seamNode)
            {
                BSIPC_ASSERT(node.is_sequence() && node.size() == 4,
                    "seaming position should be a sequence of 4-element arrays");

                SeamPoint seamPair;
                seamPair.first = Vec2(node[0].get_value<Float>(), node[1].get_value<Float>());
                seamPair.second = Vec2(node[2].get_value<Float>(), node[3].get_value<Float>());

                seamNodes.push_back(seamPair);
            }

            info.seamPoints = seamNodes;

            //BSIPC_INFO("idx: {} -- {}", info.seamSurfIndices[0], info.seamSurfIndices[1]);
            //for (const auto& seam : info.seamCoords)
            //    BSIPC_INFO("seam: ({}, {}) -- ({}, {})", seam.first[0], seam.first[1], seam.second[0], seam.second[1]);
        }
        catch (fkyaml::exception& e)
        {
            BSIPC_WARN(e.what());
            BSIPC_INFO("Error loading rest nodes position from file \"{}\". Ignoring.", path);
        }
    }

    void SolverConfig::LoadFixedMeshInfo(const fkyaml::node& root, BSIPC_OUT std::vector<FixedMeshInfo>& infos) const
    {
        auto fixedTetMeshNode = root;
        Bool existFixedTetMesh = LoadNodeFromYamlOptional(fixedTetMeshNode, root, Config::FIXED_TET_MESHES_STR);

        if (existFixedTetMesh)
        {
            BSIPC_ASSERT(fixedTetMeshNode.is_sequence(), "\"fixed_tet_meshes\" in config is not sequence")

            for (auto& mesh : fixedTetMeshNode)
            {
                BSIPC_ASSERT(mesh.is_mapping(), "subnodes of \"fixed_tet_meshes\" in config is not mapping")

                FixedMeshInfo info;
                info.type = FixedMeshType::TETRAHEDRAL;
                LoadDataFromYaml(info.meshName, mesh, Config::NAME_STR);

                this->LoadTransformations(mesh, info.transformations);
                infos.emplace_back(info);
            }
        }

        auto fixedTrigMeshNode = root;
        Bool existFixedTrigMesh = LoadNodeFromYamlOptional(fixedTrigMeshNode, root, Config::FIXED_TRIG_MESHES_STR);

        if (existFixedTrigMesh)
        {
            BSIPC_ASSERT(fixedTrigMeshNode.is_sequence(), "\"fixed_trig_meshes\" in config is not sequence")

            for (auto& mesh : fixedTrigMeshNode)
            {
                BSIPC_ASSERT(mesh.is_mapping(), "subnodes of \"fixed_trig_meshes\" in config is not mapping")

                FixedMeshInfo info;
                info.type = FixedMeshType::TRIANGLE;
                LoadDataFromYaml(info.meshName, mesh, Config::NAME_STR);

                this->LoadTransformations(mesh, info.transformations);
                infos.emplace_back(info);
            }
        }
    }

    void SolverConfig::LoadAnimatedMeshInfo(const fkyaml::node& parent, BSIPC_OUT std::optional<AnimatedMeshInfo>& info) const
    {
        auto animatedMeshNode = parent;
        Bool existAnimatedMesh = LoadNodeFromYamlOptional(animatedMeshNode, parent, Config::ANIMATED_MESHES_STR);

        if (existAnimatedMesh)
        {
            String dirName, pattern;
            UInt length;
            Float stiffness, tolerance;
            std::vector<LinearTransformation> transformations;
            LoadDataFromYaml(dirName, animatedMeshNode, Config::DIRECTORY_STR);
            LoadDataFromYaml(pattern, animatedMeshNode, Config::PATTERN_STR);
            LoadDataFromYaml(length, animatedMeshNode, Config::LENGTH_STR);
            LoadDataFromYaml<Float>(tolerance, animatedMeshNode, Config::DIFF_STR);
            LoadDataFromYaml<Float>(stiffness, animatedMeshNode, Config::STIFFNESS_STR);
            this->LoadTransformations(animatedMeshNode, transformations);

            info = AnimatedMeshInfo{ dirName, pattern, length, stiffness, tolerance, transformations };
        }
    }

    void SolverConfig::LoadExperimentLayout(const String& path)
    {
        std::ifstream ifs(path);

        try
        {
            if (!ifs)
            {
                std::string msg = "error opening config file " + path;
                throw fkyaml::exception(msg.c_str());
            }
            fkyaml::node root = fkyaml::node::deserialize(ifs);

            BSIPC_WARN("Reading configs:");
            BSIPC_INFO(ToStr(root));
            BSIPC_WARN("==================");

            // If continueOnStep is true, more will be done in Solver::InitSimulation
            LoadDataFromYamlOptional(this->continueOnStep, root, Config::CONTINUE_ON_STEP_STR);

            if (this->continueOnStep)
                LoadDataFromYamlOptional(this->startStep, root, Config::CURRENT_STEP_STR);

            // each config must explicitly specify the surfaces
            auto surfacesInfoNode = root;
            this->surfacesInfo.clear();
            LoadNodeFromYaml(surfacesInfoNode, root, Config::BS_SURFACES_STR);
            BSIPC_ASSERT(surfacesInfoNode.is_sequence(), "\"bs_surfaces\" in config is not sequence")

            // surface infos
            this->LoadSurfaceInfos(surfacesInfoNode, this->surfacesInfo);

            // fixed meshes
            this->fixedMeshInfos.clear();
            this->LoadFixedMeshInfo(root, this->fixedMeshInfos);

            // animated meshes
            this->LoadAnimatedMeshInfo(root, this->animatedMeshInfo);

            auto seamNode = root;
            this->seamInfos.clear();
            Bool existSeams = LoadNodeFromYamlOptional(seamNode, root, Config::SEAMS_STR);
            
            if (existSeams)
            {
                BSIPC_ASSERT(seamNode.is_sequence(), "\"seams\" node must be a sequence");
                for (auto& seam : seamNode)
                {
                    SeamInfo info;

                    auto surfIdxNode = seamNode;
                    LoadNodeFromYaml(surfIdxNode, seam, Config::IDX_PAIR_STR);
                    BSIPC_ASSERT(surfIdxNode.is_sequence() && surfIdxNode.size() == 2, 
                        "\"idx_pair\" node must be a sequence of length 2");

                    for (UInt i = 0; i != 2; ++i)
                    {
                        BSIPC_ASSERT(surfIdxNode[i].is_integer(), 
                            "indices in \"idx_pair\" must be integer");
                        info.seamSurfIndices[i] = surfIdxNode[i].get_value<UInt>();
                    }

                    LoadDataFromYaml(info.seamingPointsFileName, seam, Config::SEAMING_POINTS_FILE_STR);
                    String seamPath = fmt::format("seam/{}.yaml", info.seamingPointsFileName);
                    this->LoadSeamInfo(seamPath, info);

                    this->seamInfos.emplace_back(info);
                }
            }

            LoadDataFromYamlOptional(this->maxSteps, root, Config::MAX_STEPS_STR);
            LoadDataFromYamlOptional(this->maxLSIterations, root, Config::MAX_LS_ITERATIONS_STR);
            LoadDataFromYamlOptional(this->maxNewtonIterations, root, Config::MAX_NEWTON_ITERATIONS_STR);
            LoadDataFromYamlOptional(this->density, root, Config::DENSITY_STR);

            LoadDataFromYamlOptional(this->enableContact, root, Config::ENABLE_CONTACT_STR);
            LoadDataFromYamlOptional(this->contactThreshold, root, Config::CONTACT_THRESHOLD_STR);
            LoadDataFromYamlOptional(this->groundHeight, root, Config::GROUND_HEIGHT_STR);
            LoadDataFromYamlOptional(this->timestep, root, Config::TIMESTEP_STR);
            LoadDataFromYamlOptional(this->tolerance, root, Config::RESIDUE_TOLERANCE_STR);
            LoadDataFromYamlOptional(this->disableInterpatchContact, root, Config::DISABLE_INTERPATCH_CONTACT_STR);

            String quadStr;
            LoadDataFromYamlOptional(quadStr, root, Config::QUAD_SCHEME_STR);

            this->quadScheme = QuadStrToScheme(quadStr);

            /**** Available configurations for quad_scheme:
             * local                   | Place quadrature points with order per cell
             * local-edge-dense        | local, but with more quadrature points on edges/corners (3-3, otherwise 2-2)
             * local-chessboard        | local-edge-dense, but with horizontal/vertical patterns, like chessboard (2-1 or 1-2)
             * local-chessboard-dense  | alternate form of local-chessboard, with denser/sparse patterns (2-2 or 1-1)
             * global                  | Place quadrature points with order per whole element
             * global-limited          | global, with limited number of quadrature points per cell (2-2)
             * global-edge-dense       | global-limited, but with more quadrature points on edges/corners (3-3)
             * global-chessboard       | global-edge-dense, but with horizontal/vertical patterns, like chessboard (2-1 or 1-2)
             */

            LoadDataFromYamlOptional(this->quadOrder, root, Config::QUAD_ORDER_STR);
            LoadDataFromYamlOptional(this->perPatchMesh, root, Config::PER_PATCH_MESH_STR);
            LoadDataFromYamlOptional(this->screenshot, root, Config::SCREENSHOT_STR);
            LoadDataFromYamlOptional(this->writeIPCMesh, root, Config::WRITE_IPC_MESH_STR);
            LoadDataFromYamlOptional(this->writeIterMesh, root, Config::WRITE_ITER_MESH_STR);
            LoadDataFromYamlOptional(this->cacheQuadPoints, root, Config::CACHE_QUAD_POINTS_STR);
            LoadDataFromYamlOptional(this->useCachedQuadPoints, root, Config::USE_CACHED_QUAD_POINTS_STR);

            // Stiffness
            LoadDataFromYamlOptional(this->poissonRatio, root, Config::POISSON_RATIO_STR);
            LoadDataFromYamlOptional(this->thickness, root, Config::THICKNESS_STR);
            LoadDataFromYamlOptional(this->frictionMu, root, Config::FRICTION_MU_STR);
            LoadDataFromYamlOptional(this->gravConst, root, Config::GRAVITY_STR);

            if (LoadDataFromYamlOptional(this->ssYoungsModulus, root, Config::STRETCHING_YOUNGS_MODULUS_STR))
            {
                Float shearStretchRatio = 1.;
                LoadDataFromYaml(shearStretchRatio, root, Config::SHEAR_STRECTH_RATIO_STR);
                BSIPC_PEEK(shearStretchRatio);

                this->stretchStiffness = ssYoungsModulus / (2 * (1 + this->poissonRatio));
                this->shearStiffness = shearStretchRatio * this->stretchStiffness;
            }

            if (LoadDataFromYamlOptional(this->bdYoungsModulus, root, Config::BENDING_YOUNGS_MODULUS_STR))
            {
                this->bendingStiffness = this->bdYoungsModulus * this->thickness * this->thickness * this->thickness 
                    / (12 * (1 - this->poissonRatio * this->poissonRatio));
            }

            this->lambda = (this->ssYoungsModulus * this->poissonRatio) / ((1 + this->poissonRatio) * (1 - 2 * this->poissonRatio));
            this->mu = this->ssYoungsModulus / (2 * (1 + this->poissonRatio));

            BSIPC_CRITICAL("lambda: {:.2e} // mu: {:.2e}", this->lambda, this->mu);

            LoadDataFromYamlOptional(this->seamStrength, root, Config::SEAM_STRENGTH_STR);
            BSIPC_PEEK(this->seamStrength);

            Float tempContactStiffness = 0.;
            if (LoadDataFromYamlOptional(tempContactStiffness, root, Config::CONTACT_STIFFNESS_STR))
                this->contactStiffness = tempContactStiffness;
        }
        catch (fkyaml::exception& e)
        {
            BSIPC_WARN(e.what());
            BSIPC_INFO("No config specified, using default config.");
        }

        // Save the configuration for writing step logs
        try
        {
            this->SerializeConfig();
        }
        catch (fkyaml::exception& e)
        {
            BSIPC_WARN(e.what());
            BSIPC_WARN("Error serializing solver config to Yaml node.");
        }
    }

    void SolverConfig::LoadExperimentLayout()
    {
        const String configPath = "config.yaml";
        this->LoadExperimentLayout(configPath);
    }


    void SolverConfig::SerializeConfig()
    {
        fkyaml::node node(fkyaml::node_type::MAPPING);
        node[Config::CONTINUE_ON_STEP_STR] = true;
        node[Config::MAX_STEPS_STR] = this->maxSteps;
        node[Config::MAX_LS_ITERATIONS_STR] = this->maxLSIterations;
        node[Config::MAX_NEWTON_ITERATIONS_STR] = this->maxNewtonIterations;
        node[Config::DENSITY_STR] = this->density;
        node[Config::QUAD_SCHEME_STR] = QuadSchemeToStr(this->quadScheme);
        node[Config::QUAD_ORDER_STR] = this->quadOrder;
        node[Config::PER_PATCH_MESH_STR] = this->perPatchMesh;
        node[Config::SCREENSHOT_STR] = this->screenshot;
        node[Config::STRETCHING_YOUNGS_MODULUS_STR] = this->ssYoungsModulus;
        node[Config::SHEAR_STRECTH_RATIO_STR] = this->shearStiffness / this->stretchStiffness;
        node[Config::BENDING_YOUNGS_MODULUS_STR] = this->bdYoungsModulus;
        node[Config::POISSON_RATIO_STR] = this->poissonRatio;
        node[Config::THICKNESS_STR] = this->thickness;
        node[Config::FRICTION_MU_STR] = this->frictionMu;
        node[Config::CONTACT_THRESHOLD_STR] = this->contactThreshold;
        node[Config::TIMESTEP_STR] = this->timestep;
        node[Config::RESIDUE_TOLERANCE_STR] = this->tolerance;
        node[Config::CACHE_QUAD_POINTS_STR] = this->cacheQuadPoints;
        node[Config::USE_CACHED_QUAD_POINTS_STR] = this->useCachedQuadPoints;
        node[Config::GROUND_HEIGHT_STR] = this->groundHeight;
        node[Config::WRITE_IPC_MESH_STR] = this->writeIPCMesh;
        node[Config::WRITE_ITER_MESH_STR] = this->writeIterMesh;
        node[Config::DISABLE_INTERPATCH_CONTACT_STR] = this->disableInterpatchContact;

        if (this->contactStiffness.has_value())
            node[Config::CONTACT_STIFFNESS_STR] = this->contactStiffness.value();

        std::vector<fkyaml::node> surfacesData;
        for (UInt i = 0; i != this->surfacesInfo.size(); ++i)
        {
            fkyaml::node curSurface(fkyaml::node_type::MAPPING);
            const BSSurfaceInfo& info = this->surfacesInfo[i];
            curSurface[Config::IDX_STR] = info.idx;
            curSurface[Config::RESOLUTION_STR] = { { Config::X_STR, info.resolutionX }, { Config::Y_STR, info.resolutionY } };
            curSurface[Config::LINEAR_RESOLUTION_STR] = { { Config::X_STR, info.linearResolutionX }, { Config::Y_STR, info.linearResolutionY } };
            curSurface[Config::REST_SIZE_STR] = { { Config::X_STR, info.restSizeX }, { Config::Y_STR, info.restSizeY } };
            curSurface[Config::INIT_RATIO_STR] = { { Config::X_STR, info.initRatioX }, { Config::Y_STR, info.initRatioY } };
            curSurface[Config::INIT_VELOCITY_STR] = { info.initVel[0], info.initVel[1], info.initVel[2] };

            String startConfigStr = info.startConfig == StartConfig::FLAT ? Config::FLAT_STR : Config::VERTICAL_STR;
            curSurface[Config::START_CONFIG_STR] = startConfigStr;

            if (info.irregularDomain)
                curSurface[Config::REST_NODE_FILE_STR] = info.restNodeFileName;

            if (!info.fixedNodes.empty())
            {
                std::vector<fkyaml::node> fixedNodeData;
                for (UInt fixedNode : info.fixedNodes)
                    fixedNodeData.emplace_back(fixedNode);
                curSurface[Config::FIXED_NODES_STR] = fixedNodeData;
            }

            if (!info.movingNodes.empty())
            {
                std::vector<fkyaml::node> movingData;
                for (const MovingNode& movingNode : info.movingNodes)
                {
                    fkyaml::node topNode(fkyaml::node_type::MAPPING);
                    fkyaml::node indexNode = movingNode.index;
                    std::vector<fkyaml::node> velData = { movingNode.velocity[0], movingNode.velocity[1], movingNode.velocity[2] };

                    topNode[Config::MOVING_NODE_VELOCITY_STR] = velData;
                    topNode[Config::MOVING_NODE_ID_STR] = indexNode;
                    topNode[Config::MOVING_NODE_UNTIL_FRAME] = movingNode.untilFrame;
                    movingData.emplace_back(topNode);
                }
                curSurface[Config::MOVING_NODES_STR] = movingData;
            }

            surfacesData.emplace_back(curSurface);
        }
        node[Config::BS_SURFACES_STR] = surfacesData;

        std::vector<fkyaml::node> fixedTetMeshData, fixedTrigMeshData;
        for (const FixedMeshInfo& info : this->fixedMeshInfos)
        {
            fkyaml::node curMeshNode(fkyaml::node_type::MAPPING);
            curMeshNode[Config::NAME_STR] = info.meshName;

            this->SerializeTransformations(curMeshNode, info.transformations);

            if (info.type == FixedMeshType::TETRAHEDRAL)
                fixedTetMeshData.push_back(curMeshNode);
            else if (info.type == FixedMeshType::TRIANGLE)
                fixedTrigMeshData.push_back(curMeshNode);
            else
                BSIPC_WARN("[Fallthrough] Invalid fixed mesh type in config. Ignoring current mesh");
        }

        if (!fixedTetMeshData.empty())
            node[Config::FIXED_TET_MESHES_STR] = fixedTetMeshData;
        if (!fixedTrigMeshData.empty())
            node[Config::FIXED_TRIG_MESHES_STR] = fixedTrigMeshData;

        if (!this->seamInfos.empty())
        {
            std::vector<fkyaml::node> seamData;
            for (const SeamInfo& seam : this->seamInfos)
            {
                fkyaml::node curSeamNode(fkyaml::node_type::MAPPING);
                curSeamNode[Config::IDX_PAIR_STR] = { seam.seamSurfIndices[0], seam.seamSurfIndices[1] };
                curSeamNode[Config::SEAMING_POINTS_FILE_STR] = seam.seamingPointsFileName;
                seamData.emplace_back(curSeamNode);
            }

            node[Config::SEAMS_STR] = seamData;
        }

        if (this->animatedMeshInfo.has_value())
        {
            const AnimatedMeshInfo& info = this->animatedMeshInfo.value();

            fkyaml::node animatedMeshNode(fkyaml::node_type::MAPPING);
            animatedMeshNode[Config::DIRECTORY_STR] = info.dirName;
            animatedMeshNode[Config::PATTERN_STR] = "'" + info.pattern + "'";
            animatedMeshNode[Config::LENGTH_STR] = info.length;
            animatedMeshNode[Config::DIFF_STR] = info.tolerance;
            animatedMeshNode[Config::STIFFNESS_STR] = info.stiffness;
            this->SerializeTransformations(animatedMeshNode, info.transformations);

            node[Config::ANIMATED_MESHES_STR] = animatedMeshNode;
        }

        this->configRootNode = node;
        BSIPC_INFO(ToStr(this->configRootNode));
    }

    void SolverConfig::LoadTransformations(const fkyaml::node& parent, BSIPC_OUT std::vector<LinearTransformation>& transformations) const
    {
        auto transformationNode = parent;
        Bool existTransformations = LoadNodeFromYamlOptional(transformationNode, parent, Config::TRANSFORMATIONS_STR);

        if (!existTransformations)
            return;

        BSIPC_ASSERT(transformationNode.is_sequence(), "\"transformations\" in config is not sequence")

        for (auto& transformation : transformationNode)
        {
            BSIPC_ASSERT(transformation.is_mapping(), "subnodes of \"transformations\" in config is not mapping")

            String typeStr;
            LoadDataFromYaml(typeStr, transformation, Config::TRANSFORMATION_TYPE_STR);

            TransformationType type;
            if (typeStr == Config::TRANSLATION_STR)
                type = TransformationType::Translation;
            else if (typeStr == Config::ROTATION_STR)
                type = TransformationType::Rotation;
            else if (typeStr == Config::SCALING_STR)
                type = TransformationType::Scaling;
            else
            {
                BSIPC_WARN("Invalid transformation type \"{}\" in config. Ignoring current transformation", typeStr);
                continue;
            }

            auto contentNode = transformation;
            LoadNodeFromYaml(contentNode, transformation, Config::TRANSFORMATION_CONTENT_STR);
            BSIPC_ASSERT(contentNode.is_sequence() && contentNode.size() == 3, "invalid \"content\" node data");
            Vec3 content;
            for (UInt i = 0; i < 3; ++i)
            {
                BSIPC_ASSERT(contentNode[i].is_float_number(), "\"content\" node has non-floating number data");
                content[i] = contentNode[i].get_value<Float>();
            }

            transformations.emplace_back(type, content);
        }
    }

    void SolverConfig::SerializeTransformations(BSIPC_OUT fkyaml::node& parent, const std::vector<LinearTransformation>& transformations) const
    {
        if (transformations.empty())
            return;

        std::vector<fkyaml::node> transformationData;
        for (const auto& transformation : transformations)
        {
            fkyaml::node curTransformation(fkyaml::node_type::MAPPING);
            String typeStr = "";

            if (transformation.type == TransformationType::Translation)
                typeStr = Config::TRANSLATION_STR;
            else if (transformation.type == TransformationType::Rotation)
                typeStr = Config::ROTATION_STR;
            else if (transformation.type == TransformationType::Scaling)
                typeStr = Config::SCALING_STR;
            else
                BSIPC_WARN("Transformation enum fallthrough when writing to config.");

            curTransformation[Config::TRANSFORMATION_TYPE_STR] = typeStr;
            curTransformation[Config::TRANSFORMATION_CONTENT_STR] = { transformation.content[0], transformation.content[1], transformation.content[2] };

            transformationData.emplace_back(curTransformation);
        }
        parent[Config::TRANSFORMATIONS_STR] = transformationData;
    }

    std::optional<Vec3> SolverConfig::MovingNodeVelocity(UInt surfIndex, UInt nodeIndex) const
    {
        for (const MovingNode& movingNode : this->surfacesInfo[surfIndex].movingNodes)
            if (movingNode.index == nodeIndex)
                return movingNode.velocity;

        return std::nullopt;
    }
    
}
