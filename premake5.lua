workspace "B-Spline-IPC"
    architecture "x64"
    configurations {
        "Debug",   -- dev mode
        "Release"  -- release mode
    }
    startproject "B-Spline-IPC"
    
outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"

newoption {
    trigger     = "no-mkl",
    description = "Disable MKL and Pardiso usage, use Eigen's native solver instead"
}

project "B-Spline-IPC"    
    location "B-Spline-IPC"
    kind "ConsoleApp"
    language "c++"
    cppdialect "c++20"
    
    targetdir ("bin/" .. outputdir .. "/%{prj.name}")
    objdir ("bin-int/" .. outputdir .. "/%{prj.name}")

    pchheader "bspch.h"
    pchsource "%{prj.name}/src/bspch.cpp"
    
    postbuildcommands ("{COPY} dll/* ../bin/" .. outputdir .. "/%{prj.name}")

    filter "system:Windows"
        defines "BSIPC_WINDOWS"
        systemversion "latest"
        warnings "Extra"

        buildoptions {"/bigobj", "/MP"}

        staticruntime "off"                 -- set to off to avoid weird compatibility issues (/MD*)

        files {
            "%{prj.name}/src/**.h",
            "%{prj.name}/src/**.cpp",
            "%{prj.name}/ext/IPC/**.h",
            "%{prj.name}/ext/IPC/**.hpp",
            "%{prj.name}/ext/IPC/**.cpp",
            "Altrar/Altrar/src/**.h",
            "Altrar/Altrar/src/**.cpp"
        }

        excludes "Altrar/Altrar/src/EntryPoint.cpp"

        includedirs {
            "%{prj.name}/src",
            "%{prj.name}/ext",
            "Altrar/",                      -- for access from B-Spline-IPC
            "Altrar/Altrar/src",            -- for compatibility with Altrar; should be removed once Altrar switches to a static lib
            "Altrar/Altrar/src/Core",
            "Altrar/Altrar/ext",
        }

        links {
            "glfw3",
            "vulkan-1",
            "tbb",
            "amd",
            "camd",
            "colamd",
            "ccolamd",
            "cholmod",
            "openblas",
            "suitesparseconfig",
        }

        if _OPTIONS["no-mkl"] then
            libdirs {
                "Altrar/Altrar/ext/_libs",
                "%{prj.name}/lib/tbb/intel64/vc14",
                "%{prj.name}/lib/suitesparse",
                "%{prj.name}/lib/lapack",
            }

            postbuildcommands {
                "{COPY} dll/* ../bin/" .. outputdir .. "/%{prj.name}",
            }

            undefines "EIGEN_USE_MKL"
            undefines "EIGEN_USE_MKL_ALL"
        else
            MKLROOT = os.getenv("MKLROOT") or "C:/Program Files (x86)/Intel/oneAPI/mkl/latest"
            CMPLRROOT = os.getenv("CMPLR_ROOT") or "C:/Program Files (x86)/Intel/oneAPI/compiler/latest"

            includedirs {
                MKLROOT .. "/include"
            }

            libdirs {
                "Altrar/Altrar/ext/_libs",
                "%{prj.name}/lib/tbb/intel64/vc14",
                "%{prj.name}/lib/suitesparse",
                "%{prj.name}/lib/lapack",
                -- ADDED: MKL library directory
                MKLROOT .. "/lib",
                CMPLRROOT .. "/lib"
            }

            links {
                "mkl_intel_lp64",
                "mkl_intel_thread",
                "mkl_core",
                "libiomp5md"
            }

            postbuildcommands {
                "{COPY} dll/* ../bin/" .. outputdir .. "/%{prj.name}",
                "{COPY} \"" .. MKLROOT .. "/bin/*.dll\" ../bin/" .. outputdir .. "/%{prj.name}",
                "{COPY} \"" .. CMPLRROOT .. "/bin/*.dll\" ../bin/" .. outputdir .. "/%{prj.name}",
            }

            defines "EIGEN_USE_MKL_ALL"
            defines "BSIPC_USE_MKL"
        end

        disablewarnings {
            "4458"                          -- declaration of 'identifier' hides class member. Often in setters/inits this is necessary (and will not necessarily cause problems)
        }

    filter "system:Linux"
        defines "BSIPC_LINUX"
        warnings "Off"

        -- Linux simulation disables all rendering parts
        files {
            "%{prj.name}/ext/IPC/**.h",
            "%{prj.name}/ext/IPC/**.hpp",
            "%{prj.name}/ext/IPC/**.cpp",
            "%{prj.name}/src/**.h",
            "%{prj.name}/src/**.cpp",
        }

        includedirs {
            "%{prj.name}/src",
            "%{prj.name}/ext",
        }

        links {
            "tbb",
            "amd",
            "camd",
            "colamd",
            "ccolamd",
            "cholmod",
            "suitesparseconfig",
        }

        TBBROOT = os.getenv("TBBROOT")

        libdirs {
            TBBROOT,
        }

        if _OPTIONS["no-mkl"] then
            undefines "EIGEN_USE_MKL"
            undefines "EIGEN_USE_MKL_ALL"
        else
            MKLROOT = os.getenv("MKLROOT")
            CMPLRROOT = os.getenv("CMPLR_ROOT")

            includedirs {
                MKLROOT .. "/include"
            }

            libdirs {
                MKLROOT .. "/lib",
                CMPLRROOT .. "/lib",
            }

            links {
                "mkl_intel_lp64",
                "mkl_intel_thread",
                "mkl_core",
                "iomp5"
            }

            defines "BSIPC_USE_MKL"
            defines "EIGEN_USE_MKL"
            defines "EIGEN_USE_MKL_ALL"
        end

    filter "action:vs*"
        openmp "On"

    filter "files:Altrar/**"
        flags {"NoPCH"}

    filter "files:B-Spline-IPC/ext/IPC/**"
        flags {"NoPCH"}
        warnings "Off"

    filter "configurations:Debug"
        runtime "Debug"                 -- /MDd
        defines "ATR_DEBUG"
        defines "BSIPC_DEBUG"
        symbols "on"

    filter "configurations:Release"
        runtime "Release"               -- /MD
        defines "ATR_RELEASE"
        defines "BSIPC_RELEASE"
        optimize "on"
