
#ifndef FEASST_CONFINEMENT_SPHERE_H_
#define FEASST_CONFINEMENT_SPHERE_H_

#include "confinement/include/shape.h"
#include "utils/include/arguments.h"

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

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_SPHERE_H_
