#include "./barrier_planar.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

BarrierPlanar::BarrierPlanar(const double coord, const int direction,
  const int dimen)
  : coord_(coord),
    direction_(direction),
    dimen_(dimen) {
  className_.assign("BarrierPlanar");
}

double BarrierPlanar::potential(const vector<double> coordinate,
    const double diameter) {
  ASSERT(static_cast<int>(coordinate.size()) > dimen_, 
    "the requestest barrier dimension "
    << "is not present in the coordinates provided to the the potential "
    << "function");
    
  if (direction_*coordinate[dimen_] + diameter/2. > direction_*coord_) {
    return std::numeric_limits<double>::max()/1e10;
  }
  return 0;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

