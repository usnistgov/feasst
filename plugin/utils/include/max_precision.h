#ifndef FEASST_UTILS_MAX_PRECISION_H_
#define FEASST_UTILS_MAX_PRECISION_H_

#include <iomanip>  // setprecision
#include <limits>  // numeric_limits

namespace feasst {

/// Used to output maximum precision to screen
/// [e.g., INFO(MAX_PRECISMAX_PRECISIONN << "energy: " << energy)].
#define MAX_PRECISION std::setprecision(std::numeric_limits<double>::digits10+2)

}  // namespace feasst

#endif  // FEASST_UTILS_MAX_PRECISION_H_
