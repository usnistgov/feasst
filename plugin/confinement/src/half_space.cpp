#include <cmath>
#include "confinement/include/half_space.h"

namespace feasst {

class MapHalfSpace {
 public:
  MapHalfSpace() {
    HalfSpace().deserialize_map()["HalfSpace"] = MakeHalfSpace();
  }
};

static MapHalfSpace mapper_ = MapHalfSpace();

HalfSpace& HalfSpace::set_direction(
    const double direction) {
  ASSERT(std::abs(direction) > 1e-15, "direction cannot be infinitesimal");
  if (direction > 0) {
    direction_ = 1.;
  } else {
    direction_ = -1.;
  }
  return *this;
}

double HalfSpace::nearest_distance(const Position& point) const {
  DEBUG("c " << point.coord(dimension_) << " " <<
                -1.*direction_*(point.coord(dimension_) - intersection_));
  return -1.*direction_*(point.coord(dimension_) - intersection_);
}

}  // namespace feasst
