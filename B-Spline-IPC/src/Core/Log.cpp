#include "bspch.h"

#include "Log.h"

namespace BSIPC
{
	void Log::Init()
	{
		spdlog::set_pattern("%^[%T] %n: %v%$");

		Log::Logger = spdlog::stdout_color_mt("BSIPC");
		Log::Logger->set_level(spdlog::level::trace);
	}

}
