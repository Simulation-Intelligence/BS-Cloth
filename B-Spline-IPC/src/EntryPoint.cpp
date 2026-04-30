#include "bspch.h"

#include <iostream>

#include "Utility/Profiler.h"
#include "Utility/Utility.h"
#include "Utility/Exporter/ObjExporter.h"

#include "Geometry/BSpline/Primitives.h"

#include "Solver/Solver.h"
#include "Solver/BSSurfaceMgr.h"

//#define BSIPC_YAML2OBJ
//#undef BSIPC_YAML2OBJ

//#define BSIPC_DISABLE_RENDERER
//#undef BSIPC_DISABLE_RENDERER

#if defined BSIPC_YAML2OBJ
#include <filesystem>
namespace fs = std::filesystem;
#endif

int main()
{
    try
    {
        BSIPC::Log::Init();
        
        Eigen::initParallel();
#if defined BSIPC_WINDOWS
        Eigen::setNbThreads(24);            // Should be the number of physical cores of CPU
#elif defined BSIPC_LINUX
        Eigen::setNbThreads(32);            // Should be the number of physical cores of CPU
#endif

#if defined BSIPC_YAML2OBJ
        using BSIPC::UInt, BSIPC::Bool, BSIPC::String;

        BSIPC::OS::DeleteDir("objs");
        BSIPC::OS::CreateDir("objs");

        fs::path yamlStepsFolder{ "steps" };
        const UInt stepsCnt = 800;
        const String filenamePrefix = "step_";
        const UInt startIndex = 1;

        const fs::path startYamlPath = yamlStepsFolder / "start.yaml";

        BSIPC::SolverConfig config = BSIPC::SolverConfig::DefaultInfo();
        config.LoadExperimentLayout(startYamlPath.string());
        Bool sepPatch = false, renderCtrlPts = false;
        UInt surfCnt = config.surfacesInfo.size();
        BSIPC::ObjExporter exporter(sepPatch, renderCtrlPts, surfCnt);

        auto ConvertFile = [&](const fs::path path, String frameIndexStr) {
            BSIPC::BSSurfaceMgr mgr;
            mgr.ParseFromConfig(config);
            BSIPC::Solver solver(config);
            mgr.LinkSolver(solver);

            UInt surfCnt = solver.GetTargetCnt();
            std::vector<std::pair<UInt, UInt>> bsTrigSizes(surfCnt);
            for (UInt idx = 0; idx != surfCnt; ++idx)
            {
                auto surface = *(solver.GetTarget(idx));
                UInt uSize = static_cast<UInt>(surface.UDomainMax() * config.perPatchMesh + 1);
                UInt vSize = static_cast<UInt>(surface.VDomainMax() * config.perPatchMesh + 1);
                bsTrigSizes[idx] = { uSize, vSize };
            }

            solver.InitPatchOrderMesh(bsTrigSizes);

            solver.ReadFromStepLog(path.string());
            UInt meshInterval = config.perPatchMesh * config.perPatchMesh * 2;

            for (UInt idx = 0; idx != surfCnt; ++idx)
            {
                const String meshIndexPrefix = "s" + std::to_string(idx) + "_";
                const String objPrefix = meshIndexPrefix + filenamePrefix;

                exporter.Clear();
                exporter.IncludeMesh(solver.GetPatchOrderMesh(idx));
                const std::vector<BSIPC::Vec3>& ctrlPts = solver.GetControlPoints(idx);

                if (sepPatch)
                    exporter.SetMeshInterval(meshInterval);
                if (renderCtrlPts)
                    exporter.IncludeCtrlPts(ctrlPts);

                exporter.Serialize(objPrefix, idx, frameIndexStr, 0.002);
            }
        };

        ConvertFile(startYamlPath.string(), "start");

        for (UInt frameIndex = startIndex; frameIndex < startIndex + stepsCnt; ++frameIndex)
        {
            BSIPC_INFO("Loading Frame {}", frameIndex);

            const String yamlFilename = filenamePrefix + std::to_string(frameIndex) + ".yaml";
            fs::path yamlConfigPath = yamlStepsFolder / yamlFilename;
            if (fs::exists(yamlConfigPath))
                ConvertFile(yamlConfigPath, std::to_string(frameIndex));
            else
                BSIPC_ERROR("File {} does not exist", yamlConfigPath.string());
        }

        //for (UInt iterIndex = 0; iterIndex < 100; ++iterIndex)
        //{
        //    BSIPC_INFO("Loading Frame {}", iterIndex);

        //    const String yamlFilename = fmt::format("step_0_iter_{}.yaml", iterIndex);
        //    fs::path yamlConfigPath = yamlStepsFolder / yamlFilename;
        //    if (fs::exists(yamlConfigPath))
        //        ConvertFile(yamlConfigPath, std::to_string(iterIndex));
        //    else
        //        BSIPC_ERROR("File {} does not exist", yamlConfigPath.string());
        //}

#else
        BSIPC::SolverConfig config = BSIPC::SolverConfig::DefaultInfo();
        config.LoadExperimentLayout();

        BSIPC::BSSurfaceMgr mgr;
        mgr.ParseFromConfig(config);
        BSIPC::Solver solver(config);
        mgr.LinkSolver(solver);

        BSIPC::OS::DeleteFile("log.txt");

#if not defined BSIPC_DISABLE_RENDERER
        ATR::Renderer renderer;
        renderer.SetScreenshot(config.screenshot);
        renderer.ResumeAtFrame(config.startStep);

        solver.BindRenderer(&renderer);
#endif
        solver.InitSimulation();

#if not defined BSIPC_DISABLE_RENDERER
        renderer.InitRenderer();
        renderer.InitVulkan();
        if (!config.continueOnStep)
            renderer.Update();
#endif
        //try
        {
            BSIPC::Profiler::Start("siml");
            BSIPC::Profiler::Start("render");

            BSIPC::UInt stepIdx = 1;
            BSIPC::Float simlTimeSum = 0.;

            //ATR::OS::Execute("powershell -command \"Start-Sleep -Milliseconds 50\"");
#if not defined BSIPC_DISABLE_RENDERER
            while (!renderer.ShouldClose())
#else
            while (true)
#endif
            {
                BSIPC::Profiler::Update("siml");
                BSIPC::Bool stop = solver.StepForward();
                BSIPC::Profiler::Check("siml");

#if not defined BSIPC_DISABLE_RENDERER
                BSIPC::Profiler::Update("render");
                renderer.Update();
                BSIPC::Profiler::Check("render");
                BSIPC_WARN("[Step {}] Siml: {}, Render: {}", solver.CurStep(),
                    BSIPC::Profiler::GetTime("siml").value(),
                    BSIPC::Profiler::GetTime("render").value());

                //BSIPC::Float curSimlTime = BSIPC::Profiler::GetTime("siml").value();
                //simlTimeSum += curSimlTime;
                //BSIPC_WARN("Summed Simulation Time: {}", simlTimeSum);

                //std::ofstream ofs;
                //ofs.open("log.txt", std::ios::app | std::ios::out);
                //ofs << fmt::format("## Step: {}, Current Step Time: {}, Summed Time: {}\n", stepIdx, curSimlTime, simlTimeSum);
                //ofs << "========================================\n";
                //ofs.close();

                ++stepIdx;
#endif

                if (stop)
                    break;
                //ATR::OS::Execute("powershell -command \"Start-Sleep -Milliseconds 50\"");
            }
        }
        //catch (const std::exception& e)
        //{
        //    BSIPC_CRITICAL("{}", e.what());
        //}
#endif

    }
    catch (const BSIPC::Exception& e)
    {
       BSIPC_ERROR("{}", e.Info());
    }
}
