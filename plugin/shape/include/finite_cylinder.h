
#ifndef FEASST_SHAPE_FINITE_CYLINDER_H_
#define FEASST_SHAPE_FINITE_CYLINDER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "shape/include/formula_sine_wave.h"
#include "shape/include/shape_intersect.h"

namespace feasst {

/**
  While Cylinder extends infinitely, FiniteCylinder has circular end-caps,
  whose centers (along the axis of rotational symmetry) are given by two points.
 */
class FiniteCylinder : public ShapeIntersect {
 public:
  //@{
  /** @name Arguments
    - Cylinder arguments.
   */
  explicit FiniteCylinder(argtype args);
  explicit FiniteCylinder(argtype * args);
  //@}
};

inline std::shared_ptr<FiniteCylinder> MakeFiniteCylinder(argtype args) {
  return std::make_shared<FiniteCylinder>(args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_FINITE_CYLINDER_H_
