#pragma once
#include "bsfwd.h"

namespace BSIPC
{
#if defined BSIPC_WINDOWS
#undef DeleteFile
#endif
    class OS
    {
    private:
        static void Execute(String);

    public:
        static void DeleteDir(String);
        static void CreateDir(String);
        static void DeleteFile(String);
    };
}