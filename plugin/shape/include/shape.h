
#ifndef FEASST_SHAPE_SHAPE_H_
#define FEASST_SHAPE_SHAPE_H_

#include <memory>
#include <map>
#include <sstream>
#include "utils/include/arguments.h"
#include "math/include/position.h"

namespace feasst {

class Random;

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

  virtual double surface_area() const;

  virtual double volume() const;

  // HWH alternatively, rotate approximate uniform spherical grid
  /**
    Given a material of arbitrary shape that interacts with an attraction that
    scales as \f$U \approx r^{-\alpha\}f$,
    where \f$r\f$ is the radial distance from a point
    and \f$U\f$ is the interaction energy,
    integrate this interaction up to a maximum cutoff distance.

    Numerically, this is obtained by placement of uniformly random points on a
    spherical surface.
    Spherical surfaces are considered from values of \f$r\f$ from
    max_radius/num_radius to max_radius.

    Random points are used to prevent bias from approximate methods that
    attempt to generate equidistance points on sphere surfaces.

    args:
    - invert: consider the inverse shape (default: true)
    - alpha: exponential parameter (default: 6)
    - max_radius: maximum radial extent
    - num_radius: number of radius slices
    - density: number of points per unit area for each slice
   */
  double integrate(
    const Position& point,
    Random * random,
    const argtype& args = argtype()) const;

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

#endif  // FEASST_SHAPE_SHAPE_H_
