
#ifndef FEASST_SHAPE_HALF_SPACE_H_
#define FEASST_SHAPE_HALF_SPACE_H_

#include "shape/include/shape.h"
#include "utils/include/arguments.h"

namespace feasst {

/**
  A half space divides space by a plane (or line in 2D).
  If you are on the "right" side of this divding surface, then you are inside.
  Otherwise, you are not.

  The current implementation is optimized such that the dividing surface is
  assumed perpendicular to one of the coordinate axes.
  For a more arbitrarily-oriented plane, see HalfSpaceTilted.

  Thus, there are only three variables required to specify the half space.

  1. The dimension of the axis which is perpendicular to the dividing surface

  2. The coordinate value of the intersection of the dividing surface with this
     axis.

  3. The direction along this axis which is 'inside' the shape.
 */
class HalfSpace : public Shape {
 public:
  /**
    args:
    - dimension: Set the dimension of the axis which is perpendicular to the
      divider.
    - intersection: Set the value where this axis intersects the dividing
      surface.
    - direction: Set the direction at the intersection which is inside.
      The only accepted values are "1" or "-1".
   */
  explicit HalfSpace(argtype args = argtype());
  explicit HalfSpace(argtype * args);

  /// Return dimension argument.
  int dimension() const { return dimension_; }

  /// Return the intersection argument.
  double intersection() const { return intersection_; }

  /// Return dimension argument.
  int direction() const { return direction_; }

  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<HalfSpace>(istr); }
  explicit HalfSpace(std::istream& istr);
  virtual ~HalfSpace() {}

 protected:
  void serialize_half_space_(std::ostream& ostr) const;

 private:
  int dimension_;
  double intersection_;
  int direction_;
};

inline std::shared_ptr<HalfSpace> MakeHalfSpace(argtype args = argtype()) {
  return std::make_shared<HalfSpace>(args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_HALF_SPACE_H_
