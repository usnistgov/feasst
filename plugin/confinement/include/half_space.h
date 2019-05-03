
#ifndef FEASST_CONFINEMENT_HALF_SPACE_H_
#define FEASST_CONFINEMENT_HALF_SPACE_H_

#include "confinement/include/shape.h"
#include "core/include/arguments.h"

namespace feasst {

/**
  A half space divides space by a plane (or line in 2D).
  If you are on the "right" side of this divding surface, then you are inside.
  Otherwise, you are not.

  With the current implementation, the dividing surface is assumed perpendicular
  to one of the coordinate axes.

  Thus, there are only three variables required to specify the half space.

  1. The dimension of the axis which is perpendicular to the dividing surface

  2. The coordinate value of the intersection of the dividing surface with this
     axis.

  3. The direction along this axis which is 'inside' the shape.
 */
class HalfSpace : public Shape {
 public:
  HalfSpace(
    /**
      dimension : Set the dimension of the axis which is perpendicular to the
        divider.

      intersection : Set the value where this axis intersects the dividing
        surface.

      direction : Set the direction at the intersection which is inside.
        The only accepted values are "1" or "-1".
     */
    const argtype &args = argtype());

  /// Return dimension argument.
  int dimension() const { return dimension_; }

  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(376, ostr);
    feasst_serialize(dimension_, ostr);
    feasst_serialize(intersection_, ostr);
    feasst_serialize(direction_, ostr);
  }

  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<HalfSpace>(istr); }

  HalfSpace(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(376 == version, version);
    feasst_deserialize(&dimension_, istr);
    feasst_deserialize(&intersection_, istr);
    feasst_deserialize(&direction_, istr);
  }

  virtual ~HalfSpace() {}

 private:
  const std::string class_name_ = "HalfSpace";
  int dimension_;
  double intersection_;
  int direction_;
  Arguments args_;
};

inline std::shared_ptr<HalfSpace> MakeHalfSpace(
    const argtype& args = argtype()) {
  return std::make_shared<HalfSpace>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_HALF_SPACE_H_
