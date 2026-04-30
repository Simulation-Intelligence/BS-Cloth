#include "atrpch.h"

#include "ATROSSpec.h"

#include <algorithm>

namespace ATR
{

#if defined _WIN32

    namespace ATR_WINAPI
    {
// avoid macro redefinition
#undef APIENTRY

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
    }

	void OS_ChangeConsoleColor(unsigned int color)
    {
        ATR_WINAPI::HANDLE hConsole = ATR_WINAPI::GetStdHandle(static_cast<ATR_WINAPI::DWORD>(-11));            // -11 = STD_OUTPUT_HANDLE
        ATR_WINAPI::SetConsoleTextAttribute(hConsole, static_cast<ATR_WINAPI::WORD>(color));
    }

#endif

    void OS::Execute(String cmd)
    {
#if defined _WIN32
        // Since Windows requires "\" for file systems, and it needs to be escaped, first insert \\ for every "\"
        ATR_LOG_ACTION(("Executing Command: " + cmd))
        //std::replace(cmd.begin(), cmd.end(), '/', '\\');
        system(cmd.c_str());
        ATR_LOG_ACTION_END
#endif
    }

}