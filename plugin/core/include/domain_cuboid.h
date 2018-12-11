
#ifndef FEASST_CORE_DOMAIN_CUBOID_H_
#define FEASST_CORE_DOMAIN_CUBOID_H_

#include <vector>
#include "core/include/domain.h"

namespace feasst {

/**
  A cuboid-shaped domain may have unique side lengths but the angles between
  the edges of the domain are fixed at ninety degrees.
  Note that this cuboid reduces to a rectangle in 2D.
 */
class DomainCuboid : public Domain {
 public:
  /// Set the cubic box length.
  /// Return self for chain setting.
  DomainCuboid& set_cubic(const double box_length);

  // see base class for comment
  Position shift(const Position& position) const override;

  // see base class for comment
  Position random_position(Random * random) const override {
    DEBUG("side_length_ " << side_length_.str());
    return random->position_in_cuboid(side_length_);
  }

  // HWH optimization
  void wrap_opt(const Position& pos1, const Position& pos2, Position * rel, double * r2) const {
    *r2 = 0;
    const std::vector<double>& side = side_length_.coord();
    std::vector<double>* dxv = (*rel).get_coord();
    for (int dim = 0; dim < 3; ++dim) {
      (*dxv)[dim] = pos1.coord()[dim] - pos2.coord()[dim];
      if ((*dxv)[dim] >  0.5*side[dim]) {
        (*dxv)[dim] -= side[dim];
      } else if ((*dxv)[dim] < -0.5*side[dim]) {
        (*dxv)[dim] += side[dim];
      }
      *r2 += (*dxv)[dim]*(*dxv)[dim];
    }
  }

  virtual ~DomainCuboid() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_DOMAIN_CUBOID_H_
