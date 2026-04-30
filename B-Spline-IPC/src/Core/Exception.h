#pragma once

#include "Type.h"

namespace BSIPC
{
    class Exception
    {
    public:
        Exception(const String& msg, const String& file, const Int& line, const String& function) : 
            msg(msg), file(file), line(line), function(function) {  }

        inline const String What() const { return this->msg; }
        inline const String Info() const
        {
            String info = "Exception: ";
            info += this->msg;

            String debugMsg = "\nFile: " + this->file + "\nLine: " + std::to_string(this->line) + ", Function: " + this->function;

            return info + debugMsg;
        }

    private:
        String msg;
        String file;
        Int line;
        String function;
    };
}

#define BSIPC_ENABLE_ASSERT

#if defined BSIPC_ENABLE_ASSERT
#define BSIPC_ASSERT(x, msg) \
    if (!(x)) \
        throw ::BSIPC::Exception(msg, __FILE__, __LINE__, __FUNCTION__);
#else
#define BSIPC_ASSERT(x, msg)
#endif
