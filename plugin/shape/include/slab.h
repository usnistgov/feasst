
#ifndef FEASST_SHAPE_SLAB_H_
#define FEASST_SHAPE_SLAB_H_

#include <memory>
#include "utils/include/arguments.h"
#include "shape/include/shape_intersect.h"

namespace feasst {

/**
  A slab is the intersection of two half spaces.
 */
class Slab : public ShapeIntersect {
 public:
  /**
    args:
    - dimension: The slab surface is perpendicular to this dimensional axis.
    - bound0: Set a lower or upper value of the slab.
    - bound1: Set the second bound, upper or lower, respectively.
   */
  explicit Slab(argtype args = argtype());
};

inline std::shared_ptr<Slab> MakeSlab(argtype args = argtype()) {
  return std::make_shared<Slab>(args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_SLAB_H_
