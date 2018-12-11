
#include <vector>
#include "core/include/domain_cuboid.h"

namespace feasst {

DomainCuboid& DomainCuboid::set_cubic(const double box_length) {
  std::vector<double> cube = {box_length, box_length, box_length};
  Position side_length;
  side_length.set_vector(cube);
  set_side_length(side_length);
  return *this;
}

Position DomainCuboid::shift(const Position& position) const {
  Position displacement(position);
  for (int dim = 0; dim < position.size(); ++dim) {
    if (periodic_[dim]) {
      const double box_length = side_length().coord(dim);
      const double half_box_length = 0.5*box_length;
      const double coord = position.coord(dim);
      if (coord > half_box_length) {
        displacement.set_coord(dim, -box_length);
      } else if (coord < -half_box_length) {
        displacement.set_coord(dim, box_length);
      } else {
        displacement.set_coord(dim, 0);
      }
    } else {
      displacement.set_coord(dim, 0);
    }
  }
  return displacement;
}

}  // namespace feasst
