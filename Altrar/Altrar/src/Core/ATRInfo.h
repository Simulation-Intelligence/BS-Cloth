#pragma once

// Defines the logger help functions in ATR
// All Logging functions should be inlined as this header is included multiple times throughout the project

#include "ATRType.h"
#include "ATROSSpec.h"

#include <iostream>
#include <iomanip>

namespace ATR
{

#if not defined ATR_OUTPUT_SILENT

#define ATR_FILL_LEFT 20
#define ATR_FILL_LENGTH 80

#define ATR_FILL_STR_WITH(fill, leftsize, size, content)\
	std::left << std::setfill(fill) << std::setw(leftsize) << "" << std::left << std::setfill(fill) << std::setw(size) << content

#define ATR_PRINT(x) {std::cout << x << std::endl;}

// TODO change this to spdlog
#define ATR_LOG(x) {std::cout << "[Log] >> " << x << std::endl;}
#define ATR_LOG_SUB(x) {std::cout << "[Log] >> " << Format::item << x << std::endl;}
#define ATR_LOG_PART(name) {std::cout << ATR_FILL_STR_WITH('=', ATR_FILL_LEFT, ATR_FILL_LENGTH, name) << std::endl;}
#define ATR_LOG_SECTION(name) {std::cout << ATR_FILL_STR_WITH('-', ATR_FILL_LEFT, ATR_FILL_LENGTH, name) << std::endl;} 
#define ATR_LOG_ACTION(name) {std::cout << ATR_FILL_STR_WITH('#', 5, ATR_FILL_LENGTH, " " + name + "  ") << std::endl;}
#define ATR_LOG_ACTION_END {std::cout << ATR_FILL_STR_WITH('#', ATR_FILL_LEFT, ATR_FILL_LENGTH, "") << std::endl;}

#if defined ATR_VERBOSE
	#define ATR_LOG_VERBOSE(x) {std::cout << "[Log] >> " << x << std::endl;}
	#define ATR_PRINT_VERBOSE(x) {std::cout << x << std::endl;}
#else
    #define ATR_LOG_VERBOSE(x) {}
	#define ATR_PRINT_VERBOSE(x) {}
#endif

#define ATR_ERROR(x) { \
	OS_ChangeConsoleColor(12); \
	std::cout << "[ERROR] "; \
	OS_ChangeConsoleColor(7); \
	std::cout << x << std::endl;}

// The Logged `x` must be of type which has operator<< defined. Do not append the lagging '\n' to the output of `x`
#define ATR_PEEK(x) {std::cout << "| Variable: " << #x << " " << "| Type: " << typeid(x).name() << std::endl;\
	std::cout << "| Value:" << std::endl << x << std::endl; }

#endif

	struct Format
	{
		static inline String subitem = "  - ";
		static inline String item = "- ";
	};

}
