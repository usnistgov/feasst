
#ifndef FEASST_SHAPE_CYLINDER_H_
#define FEASST_SHAPE_CYLINDER_H_

#include "shape/include/shape.h"
#include "utils/include/arguments.h"

namespace feasst {

/**
  A cylinder is given by an axis of rotational symmetry and a radius.
  The axis is described by two points.
 */
class Cylinder : public Shape {
 public:
  /**
    args:
    - radius: radius of the cylinder.
    - first_point: set the unique key for the first_point positions.
      Thus, arguments of "key[i]" are expected to follow.
      The "[i]" is to be substituted for integer dimensions 0, 1, 2, ...
      The "[i]" are also expected to be in order, starting from 0.
    - second_point: as described for first_point.
   */
  explicit Cylinder(argtype args);
  explicit Cylinder(argtype * args);

  const Position& first_point() const { return point0_; }
  const Position& second_point() const { return point1_; }

  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<Cylinder>(istr); }
  std::shared_ptr<Shape> create(argtype * args) const override {
    return std::make_shared<Cylinder>(args); }
  explicit Cylinder(std::istream& istr);
  virtual ~Cylinder() {}

 private:
  double radius_;
  Position point0_, point1_;
};

inline std::shared_ptr<Cylinder> MakeCylinder(argtype args) {
  return std::make_shared<Cylinder>(args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_CYLINDER_H_
