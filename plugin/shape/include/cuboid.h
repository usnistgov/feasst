
#ifndef FEASST_SHAPE_CUBOID_H_
#define FEASST_SHAPE_CUBOID_H_

#include "shape/include/shape.h"
#include "utils/include/arguments.h"

namespace feasst {

class Random;

/**
  A cuboid is given by a center point and side lengths in each dimension.
  Not fully implemented.
  Implement as the intersection of 3 perpendicular slabs
 */
class Cuboid : public Shape {
 public:
  /// Construct given side lengths and center.
  Cuboid(const Position& side_lengths, const Position& center);

  /// Same as above, but the center is assumed to be the origin (in 3D).
  Cuboid(const Position& side_lengths)
    : Cuboid(side_lengths, Position({0, 0, 0})) {}

  /// Same as above, but all the side lengths are the same (cube, in 3D).
  Cuboid(const double cubic_side_length)
    : Cuboid(Position({cubic_side_length, cubic_side_length, cubic_side_length})) {}

  double nearest_distance(const Position& point) const override;
  double surface_area() const override;
  double volume() const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<Cuboid>(istr); }
  explicit Cuboid(std::istream& istr);
  virtual ~Cuboid() {}

 private:
  Position side_lengths_, center_;
};

inline std::shared_ptr<Cuboid> MakeCuboid(
    const Position& side_lengths,
    const Position& center) {
  return std::make_shared<Cuboid>(side_lengths, center);
}

inline std::shared_ptr<Cuboid> MakeCuboid(const Position& side_lengths) {
  return std::make_shared<Cuboid>(side_lengths);
}

inline std::shared_ptr<Cuboid> MakeCuboid(const double cubic_side_length) {
  return std::make_shared<Cuboid>(cubic_side_length);
}

inline std::shared_ptr<Cuboid> MakeCube(const double side_length) {
  return MakeCuboid(side_length);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_CUBOID_H_
