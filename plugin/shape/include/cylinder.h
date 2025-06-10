
#ifndef FEASST_SHAPE_CYLINDER_H_
#define FEASST_SHAPE_CYLINDER_H_

#include "shape/include/shape.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  A cylinder is given by an axis of rotational symmetry and a radius.
  The axis is described by two points.
 */
class Cylinder : public Shape {
 public:
  //@{
  /** @name Arguments
    - radius: radius of the cylinder.
    - first_point: comma-separated values for the positions in each dimension.
    - second_point: as described for first_point.
   */
  explicit Cylinder(argtype args);
  explicit Cylinder(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

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

  //@}
 private:
  double radius_;
  Position point0_, point1_;
};

inline std::shared_ptr<Cylinder> MakeCylinder(argtype args) {
  return std::make_shared<Cylinder>(args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_CYLINDER_H_
