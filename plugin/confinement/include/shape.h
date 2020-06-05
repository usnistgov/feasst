
#ifndef FEASST_CONFINEMENT_SHAPE_H_
#define FEASST_CONFINEMENT_SHAPE_H_

#include <memory>
#include <map>
#include <sstream>
#include "math/include/position.h"

namespace feasst {

/**
  Shapes may be defined by either a simple mathematical
  formula or an interpolated data table.
 */
class Shape {
 public:
  /// Return the distance from the point to the nearest point on the surface.
  /// The distance is negative if the point is inside of the shape and positive
  /// if it is outside.
  virtual double nearest_distance(const Position& point) const = 0;

  /// Return true if the point is inside of the shape.
  bool is_inside(const Position& point) const;

  /// Return true if the sphere of given center point and diameter is entirely
  /// inside of the shape.
  bool is_inside(const Position& point, const double diameter) const;

  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Shape> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Shape> >& deserialize_map();
  std::shared_ptr<Shape> deserialize(std::istream& istr);
  virtual ~Shape() {}
};

// An object which contains a shape.
class ShapedEntity {
 public:
  ShapedEntity() {}
  ShapedEntity(std::shared_ptr<Shape> shape) { shape_ = shape; }

  /// Return the shape.
  const std::shared_ptr<Shape> shape() const { return shape_; }

  void serialize(std::ostream& ostr) const;
  explicit ShapedEntity(std::istream& istr);

 private:
  std::shared_ptr<Shape> shape_;
};

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_SHAPE_H_
