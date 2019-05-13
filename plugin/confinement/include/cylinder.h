
#ifndef FEASST_CONFINEMENT_CYLINDER_H_
#define FEASST_CONFINEMENT_CYLINDER_H_

#include "confinement/include/shape.h"
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

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(629, ostr);
    feasst_serialize(radius_, ostr);
    feasst_serialize_fstobj(point0_, ostr);
    feasst_serialize_fstobj(point1_, ostr);
  }

  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<Cylinder>(istr); }

  Cylinder(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(629 == version, version);
    feasst_deserialize(&radius_, istr);
    feasst_deserialize_fstobj(&point0_, istr);
    feasst_deserialize_fstobj(&point1_, istr);
  }

  virtual ~Cylinder() {}

 private:
  const std::string class_name_ = "Cylinder";
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

#endif  // FEASST_CONFINEMENT_CYLINDER_H_
