
#ifndef FEASST_SHAPE_SPHERE_H_
#define FEASST_SHAPE_SPHERE_H_

#include "shape/include/shape.h"
#include "utils/include/arguments.h"

namespace feasst {

class Random;

/**
  A sphere is given by a center point and a radius.
 */
class Sphere : public Shape {
 public:
  Sphere(
    /**
      radius : Set the radius of the sphere.
     */
    const argtype &args,
    /// position of the center of the sphere.
    const Position center);

  /// Same as above, but the center is assumed to be the origin.
  Sphere(const argtype &args) : Sphere(args, Position({0, 0, 0})) {}

  double nearest_distance(const Position& point) const override;
  double surface_area() const override;
  double volume() const override;

  /**
    Generate approximately equi-distance points on the surface using golden
    spiral methodology.

    http://www.softimageblog.com/archives/115
    https://en.wikipedia.org/wiki/Tammes_problem
   */
  void surface_mesh(
    /// Give the target density of points per surface area.
    /// The actualy density may deviate slightly, and can be double checked
    /// based on the number of returned points.
    const double target_density,
    std::vector<Position> * points) const;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<Sphere>(istr); }
  explicit Sphere(std::istream& istr);
  virtual ~Sphere() {}

 private:
  const std::string class_name_ = "Sphere";
  double radius_;
  Position center_;
  Arguments args_;
};

inline std::shared_ptr<Sphere> MakeSphere(
    const argtype& args,
    const Position center) {
  return std::make_shared<Sphere>(args, center);
}

inline std::shared_ptr<Sphere> MakeSphere(const argtype& args) {
  return std::make_shared<Sphere>(args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_SPHERE_H_
