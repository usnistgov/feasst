
#ifndef FEASST_CONFINEMENT_SPHERE_H_
#define FEASST_CONFINEMENT_SPHERE_H_

#include "confinement/include/shape.h"
#include "core/include/arguments.h"

namespace feasst {

/**
  A sphere is given by a center point and a radius.
 */
class Sphere : public Shape {
 public:
  Sphere(
    /**
      radius : Set the radius of the cylinder.
     */
    const argtype &args,
    /// position of the center of the sphere.
    const Position center);

  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(629, ostr);
    feasst_serialize(radius_, ostr);
    feasst_serialize_fstobj(center_, ostr);
  }

  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<Sphere>(istr); }

  Sphere(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(629 == version, version);
    feasst_deserialize(&radius_, istr);
    feasst_deserialize_fstobj(&center_, istr);
  }

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

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_SPHERE_H_
