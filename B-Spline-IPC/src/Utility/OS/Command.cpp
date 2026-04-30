#include "bspch.h"

#include "Command.h"

#include <algorithm>
#include <iostream>

namespace BSIPC
{
    void OS::Execute(String cmd)
    {
        std::cout << ">> Executing Command: " << cmd << "\n";
        system(cmd.c_str());
    }

#if defined BSIPC_WINDOWS
    void OS::DeleteDir(String dir)
    {
        OS::Execute("rmdir /s /q " + dir);
    }

    void OS::CreateDir(String dir)
    {
        OS::Execute("mkdir " + dir);
    }

    void OS::DeleteFile(String file)
    {
        OS::Execute("del " + file);
    }

#elif defined BSIPC_LINUX

    void OS::DeleteDir(String dir)
    {
        OS::Execute("rm -r " + dir);
    }

    void OS::CreateDir(String dir)
    {
        OS::Execute("mkdir " + dir);
    }

    void OS::DeleteFile(String file)
    {
        OS::Execute("rm " + file);
    }

#endif

}