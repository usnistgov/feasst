
#ifndef FEASST_CONFINEMENT_SHAPE_H_
#define FEASST_CONFINEMENT_SHAPE_H_

#include <memory>
#include "core/include/position.h"

namespace feasst {

/**
  Shapes may be defined by either a simple mathematical
  formula or an interpolated data table.
 */
class Shape {
 public:
  /// Return the distance from the point to the nearest point on the surface.
  /// The distance is negative if the point is inside of the shape.
  virtual double nearest_distance(const Position& point) const = 0;

  /// Return true if the point is inside of the shape.
  bool is_inside(const Position& point) const;

  /// Return true if the sphere of given center point and diameter is entirely
  /// inside of the shape.
  bool is_inside(const Position& point, const double diameter) const;
};

// An object which contains a shape.
class ShapedEntity {
 public:
  /// Set the shape.
  void set_shape(const std::shared_ptr<Shape> shape) { shape_ = shape; }

  /// Return the shape.
  const std::shared_ptr<Shape> shape() const { return shape_; }

 private:
  std::shared_ptr<Shape> shape_;
};

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
  /// Set the dimension of the axis which is perpendicular to the divider.
  HalfSpace& set_dimension(const int dimension) {
    dimension_ = dimension;
    return *this;
  }

  /// Set the value where this axis intersects the dividing surface.
  HalfSpace& set_intersection(const int intersection) {
    intersection_ = intersection;
    return *this;
  }

  /// Set the direction at the intersection which is inside (e.g., positive
  /// or negative).
  HalfSpace& set_direction(const double direction);

  double nearest_distance(const Position& point) const override;

 private:
  int dimension_;
  double intersection_;
  int direction_;
};

class ShapeIntersect : public Shape {
 public:
  ShapeIntersect(const std::shared_ptr<Shape> shape1,
                 const std::shared_ptr<Shape> shape2);

  double nearest_distance(const Position& point) const override;

 private:
  const std::shared_ptr<Shape> shape1_, shape2_;
};

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_SHAPE_H_
