
#ifndef FEASST_SHAPE_SHAPE_H_
#define FEASST_SHAPE_SHAPE_H_

#include <memory>
#include <map>
#include <sstream>
#include "utils/include/arguments.h"
#include "math/include/position.h"

namespace feasst {

class Random;
class Sphere;

/**
  Shapes may be defined by either a simple mathematical
  formula or an interpolated data table.
 */
class Shape {
 public:
  Shape() {}

  /// Return the distance from the point to the nearest point on the surface.
  /// The distance is negative if the point is inside of the shape and positive
  /// if it is outside.
  virtual double nearest_distance(const Position& point) const = 0;

  /// Return true if the point is inside of the shape.
  virtual bool is_inside(const Position& point) const;

  /// Return true if the sphere of given center point and diameter is entirely
  /// inside of the shape.
  virtual bool is_inside(const Position& point, const double diameter) const;

  /// Return the surface area (if implemented).
  virtual double surface_area() const;

  /// Return the volume (if implemented).
  virtual double volume() const;

  /**
    Return a grid of points inside the surface for visualization.
    First, construct a grid of points inside a cuboid with furthest corners
    given by the upper and lower points.
    Remove points not inside the shape.
   */
  std::vector<Position> grid(const Position& upper, const Position& lower,
    /// use this many points to span each dimension
    const int num_in_each_dimen = 1e2);

  /**
    Given a material of arbitrary shape that interacts with an attraction of the
    form \f$U = \sum_i \frac{\epsilon_i}{r^{\alpha_i}}\f$,
    where \f$r\f$ is the radial distance from a point
    and \f$\epsilon\f$ is the interaction energy constant factor,
    integrate this interaction up to a maximum cutoff distance.

    Spherical shells are numerically integrated with values of \f$r\f$ from
    0 to max_radius.

    An approximate method is obtained to grid the spherical shell with
    equidistance points.
    To avoid bias, over many averages, this shell is randomly rotated.

    args:
    - invert: consider the inverse shape (default: true)
    - alpha[i]: add the i-th exponential parameter (default: 6).
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only only alpha, the "[i]" is optional.
    - epsilon[i]: add the i-th constant factor (default: -1).
      The "[i]" is as described above, and each alpha must have
      a corresponding epsilon.
    - max_radius: maximum radial extent of spherical shells.
    - num_shells: number of spherical shells.
    - points_per_slice: number of points in each spherical shell.
   */
  double integrate(
    const Position& point,
    Random * random,
    argtype args = argtype());

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Shape> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Shape> >& deserialize_map();
  std::shared_ptr<Shape> deserialize(std::istream& istr);
  virtual ~Shape() {}

 protected:
  std::string class_name_ = "Shape";
  void serialize_shape_(std::ostream& ostr) const;
  explicit Shape(std::istream& istr);

 private:
  // temporary cache
  std::shared_ptr<std::vector<std::vector<Position> > > meshes_;
  std::shared_ptr<std::vector<std::shared_ptr<Sphere> > > spheres_;
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
