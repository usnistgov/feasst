
#include <vector>
#include "core/include/domain_triclinic.h"
#include "core/include/debug.h"

namespace feasst {

Position DomainTriclinic::shift(const Position& position) const {
  ASSERT(cells().size() == 0, "cell list not implemented here");
  ASSERT((position.size() == 2) || (position.size() == 3),
    "triclinic implemented for only 2 or 3 dimensions.");
  Position displacement(position);
  displacement.set_to_origin();
  double shiftx = 0., shifty = 0., shiftz = 0.;
  if (position.size() == 3) {
    const double dz = position.coord(2);
    const double lz = side_length_.coord(2);
    if (periodic_[2] && (dz > 0.5*lz)) {
      if (dz < 0.) {
        shiftz = lz;
        shifty = yz_;
        shiftx = xz_;
      } else {
        shiftz = -lz;
        shifty = -yz_;
        shiftx = -xz_;
      }
    }
    displacement.set_coord(2, shiftz);
  }
  const double dy = position.coord(1);
  const double dx = position.coord(0);
  const double ly = side_length_.coord(1);
  const double lx = side_length_.coord(0);
  if (periodic_[1] && (dy > 0.5*ly)) {
    if (dy < 0.) {
      shifty = ly;
      shiftx = xy_;
    } else {
      shifty = -ly;
      shiftx = -xy_;
    }
  }
  if (periodic_[0] && (dx > 0.5*lx)) {
    if (dx < 0.) {
      shiftx = lx;
    } else {
      shiftx = -lx;
    }
  }
  displacement.set_coord(1, shifty);
  displacement.set_coord(0, shiftx);
  return displacement;
}

Position DomainTriclinic::random_position(Random * random) const {
  ERROR("not implemented");
  return Position();
}

}  // namespace feasst
