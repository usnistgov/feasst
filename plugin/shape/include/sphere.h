
#ifndef FEASST_SHAPE_SPHERE_H_
#define FEASST_SHAPE_SPHERE_H_

#include "shape/include/shape.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Random;

/**
  A sphere is given by a center point and a radius.
 */
class Sphere : public Shape {
 public:
  //@{
  /** @name Arguments
    - radius: Set the radius of the sphere (default: 1).
    - center: set the unique key for the center positions.
      Thus, arguments of "key[i]" are expected to follow.
      The "[i]" is to be substituted for integer dimensions 0, 1, 2, ...
      The "[i]" are also expected to be in order, starting from 0.
      If center arg is not used, a three dimensional origin is assumed.
   */
  explicit Sphere(argtype args = argtype());
  explicit Sphere(argtype * args);
  //@}
  /** @name Public Functions
   */
  //@{
  const Position& center() const { return center_; }
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
    /// The number of points.
    const int num,
    std::vector<Position> * points) const;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<Sphere>(istr); }
  std::shared_ptr<Shape> create(argtype * args) const override {
    return std::make_shared<Sphere>(args); }
  explicit Sphere(std::istream& istr);
  virtual ~Sphere() {}

  //@}
 private:
  double radius_;
  Position center_;
};

inline std::shared_ptr<Sphere> MakeSphere(argtype args = argtype()) {
  return std::make_shared<Sphere>(args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_SPHERE_H_
