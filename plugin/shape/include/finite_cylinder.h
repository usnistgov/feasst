
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
  FiniteCylinder(
    /**
      radius : Set the radius of the cylinder.
     */
    const argtype& args,
    /// the first end-point along the cylinder's axis of symmetry.
    const Position& point0,
    /// the second end-point along the cylinder's axis of symmetry.
    const Position& point1);
};

inline std::shared_ptr<FiniteCylinder> MakeFiniteCylinder(
    const argtype &args,
    const Position& point0,
    const Position& point1) {
  return std::make_shared<FiniteCylinder>(args, point0, point1);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_FINITE_CYLINDER_H_
