
#ifndef FEASST_CORE_CONSTANTS_H_
#define FEASST_CORE_CONSTANTS_H_

#include <limits>

namespace feasst {

/// An almost infinity large number to the limits of double precision
constexpr double NEAR_INFINITY = std::numeric_limits<double>::max()/1e10;

/// The smallest number within numerical precision
constexpr double NEAR_ZERO = 1e-15;

/// PI = 3.14159265358979323846264338327950 truncated to double precision
constexpr double PI = 3.14159265358979323846264338327950;

}  // namespace feasst

#endif  // FEASST_CORE_CONSTANTS_H_
