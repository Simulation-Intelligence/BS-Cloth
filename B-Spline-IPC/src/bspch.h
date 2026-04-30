#pragma once

#include "bsfwd.h"

/** Nameing Convention
 * - Classes/Function: PascalCase
 * - Variables: camelCase (except for single-letter physical variables, in consistency with the doc)
 * - Enumerations/Macros: ALL_CAPS
 *
 * Functions/Variables related to partial derivatives should contain special prefix:
 * - Functions: PD_*
 * - Variables: pd_*
 */

// BSIPC Core
#include "Core/Log.h"

// BSIPC Math Basics
#include "Geometry/BSpline/BSBasis.h"
#include "Utility/Constants.h"

// STL Tools
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <ranges>

// STL Types
#include <fstream>

// tbb
#include "oneapi/tbb.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/concurrent_vector.h"
#include "oneapi/tbb/concurrent_unordered_map.h"

// OS Utilities
#include "Utility/OS/Command.h"

#if defined BSIPC_DEBUG
#undef NDEBUG
#elif defined BSIPC_RELEASE
#define NDEBUG
#endif

#if defined BSIPC_LINUX
#define BSIPC_DISABLE_RENDERER
#endif
