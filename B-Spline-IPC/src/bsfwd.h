#pragma once

// Forward declarations: all includes used in declarations should be placed here

// All BS* files/classes are assumed to be quadratic (at least for now)

// STL Containers
#include <array>
#include <set>
#include <vector>
#include <unordered_set>

// STL Utilities
#include <optional>

// BSIPC Core
#include "Core/Type.h"              // this includes Eigen library
#include "Core/Exception.h"

// Control whether numerical derivatives are printed to log
#define BSIPC_NUMERIC_TEST
//#undef BSIPC_NUMERIC_TEST

#if defined BSIPC_LINUX
#undef BSIPC_NUMERIC_TEST           // No numeric test on experiment machines
#endif

// Control whether use linear primitives for contact
#define BSIPC_LINEAR_CONTACT

#define BSIPC_OUT
