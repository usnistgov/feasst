
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
  Cylinder(
    /**
      radius : Set the radius of the cylinder.
     */
    const argtype &args,
    /// one point on the cylinder's axis of symmetry
    const Position point0,
    /// a second point on the cylinder's axis of symmetry
    const Position point1);

  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<Cylinder>(istr); }
  explicit Cylinder(std::istream& istr);
  virtual ~Cylinder() {}

 private:
  double radius_;
  Position point0_, point1_;
  Arguments args_;
};

inline std::shared_ptr<Cylinder> MakeCylinder(
    const argtype& args,
    const Position point0,
    const Position point1) {
  return std::make_shared<Cylinder>(args, point0, point1);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_CYLINDER_H_
