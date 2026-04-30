#pragma once

#pragma warning(push, 0)
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#pragma warning(pop)

namespace BSIPC
{
	class Log
	{
	public:
		static void Init();

		static std::shared_ptr<spdlog::logger>& GetLogger() { return Log::Logger; }
	private:
		inline static std::shared_ptr<spdlog::logger> Logger;
	};
}

#define BSIPC_PERF_INFO(...)		::BSIPC::Log::GetLogger()->warn(__VA_ARGS__)

#if not defined BSIPC_PERF_TEST

#define BSIPC_TRACE(...)		    ::BSIPC::Log::GetLogger()->trace(__VA_ARGS__)
#define BSIPC_INFO(...)				::BSIPC::Log::GetLogger()->info(__VA_ARGS__)
#define BSIPC_WARN(...)				::BSIPC::Log::GetLogger()->warn(__VA_ARGS__)
#define BSIPC_ERROR(...)		    ::BSIPC::Log::GetLogger()->error(__VA_ARGS__)
#define BSIPC_CRITICAL(...)			::BSIPC::Log::GetLogger()->critical(__VA_ARGS__)
#define BSIPC_FATAL(...)		    ::BSIPC::Log::GetLogger()->fatal(__VA_ARGS__)

#define BSIPC_PEEK(var)				BSIPC_INFO("{}: {} with type [{}]", #var, var, typeid(var).name())
#define BSIPC_PEEK_EIGEN(var)		BSIPC_INFO("{}: {} with type [{}]", #var, BSIPC::ToStr(var), typeid(var).name())

#define BSIPC_INFO_SEP               BSIPC_INFO("----------------------------------------");

#endif
