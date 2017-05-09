#include "./barrier_planar.h"

namespace feasst {

BarrierPlanar::BarrierPlanar(const double coord, const int direction,
  const int dimen)
  : coord_(coord),
    direction_(direction),
    dimen_(dimen) {
  className_.assign("BarrierPlanar");
}

double BarrierPlanar::potential(const vector<double> coordinate,
    const double diameter) {
  ASSERT(coordinate.size() > dimen_, "the requestest barrier dimension "
    << "is not present in the coordinates provided to the the potential "
    << "function");
    
  if (direction_*coordinate[dimen_] + diameter/2. > direction_*coord_) {
    return std::numeric_limits<double>::max()/1e10;
  }
  return 0;
}

}  // namespace feasst

