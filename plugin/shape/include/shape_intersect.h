
#ifndef FEASST_SHAPE_SHAPE_INTERSECT_H_
#define FEASST_SHAPE_SHAPE_INTERSECT_H_

#include <memory>
#include <sstream>
#include "shape/include/shape.h"

namespace feasst {

/**
  Represents the intersection of two shapes.
  An intersection of two shapes is the region which is inside both shapes.
  This is implemented by returning the largest value of the nearest distance.
 */
class ShapeIntersect : public Shape {
 public:
  // This constructor only to be used for serialization.
  ShapeIntersect() { class_name_ = "ShapeIntersect"; }

  /// Constructor
  ShapeIntersect(std::shared_ptr<Shape> shape1,
                 std::shared_ptr<Shape> shape2);

  // Set
  void set(std::shared_ptr<Shape> shape1,
           std::shared_ptr<Shape> shape2);

  bool is_inside(const Position& point) const override;
  bool is_inside(const Position& point, const double diameter) const override;
  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<ShapeIntersect>(istr); }
  explicit ShapeIntersect(std::istream& istr);
  virtual ~ShapeIntersect() {}

 protected:
  void serialize_shape_intersect_(std::ostream& ostr) const;

 private:
  std::shared_ptr<Shape> shape1_, shape2_;
};

inline std::shared_ptr<ShapeIntersect> MakeShapeIntersect(
    std::shared_ptr<Shape> shape1,
    std::shared_ptr<Shape> shape2) {
  return std::make_shared<ShapeIntersect>(shape1, shape2);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_SHAPE_INTERSECT_H_
