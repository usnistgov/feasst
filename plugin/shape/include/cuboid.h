
#ifndef FEASST_SHAPE_CUBOID_H_
#define FEASST_SHAPE_CUBOID_H_

#include "shape/include/shape.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Random;

/**
  A cuboid is given by a center point and side lengths in each dimension.
  Not fully implemented.
  Implement as the intersection of 3 perpendicular slabs
 */
class Cuboid : public Shape {
 public:
  //@{
  /** @name Arguments
    - cubic_side_length: side length of cube.
    - side_lengths: comma-separated values for the length in each dimension.
      Cannot be used in conjunction with cubic_side_length.
    - center: comma-separated values for the positions in each dimension.
      (default: origin).
   */
  explicit Cuboid(argtype args);
  explicit Cuboid(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  double nearest_distance(const Position& point) const override;
  double surface_area() const override;
  double volume() const override;
  const Position& center() const { return center_; }

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<Cuboid>(istr); }
  std::shared_ptr<Shape> create(argtype * args) const override {
    return std::make_shared<Cuboid>(args); }
  explicit Cuboid(std::istream& istr);
  virtual ~Cuboid() {}

  //@}
 private:
  Position side_lengths_, center_;
};

inline std::shared_ptr<Cuboid> MakeCuboid(argtype args) {
  return std::make_shared<Cuboid>(args); }

}  // namespace feasst

#endif  // FEASST_SHAPE_CUBOID_H_
