
#ifndef FEASST_SHAPE_SHAPE_UNION_H_
#define FEASST_SHAPE_SHAPE_UNION_H_

#include <memory>
#include <sstream>
#include "shape/include/shape.h"

namespace feasst {

/**
  Represents the union of two shapes.
  The union of two shapes is the region that is inside either shape.
  This is done by returning the smallest value of the nearest distance.
 */
class ShapeUnion : public Shape {
 public:
  // This constructor only to be used for serialization.
  ShapeUnion() {}

  ShapeUnion(std::shared_ptr<Shape> shape1,
             std::shared_ptr<Shape> shape2);

  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<ShapeUnion>(istr); }
  explicit ShapeUnion(std::istream& istr);
  virtual ~ShapeUnion() {}

 private:
  const std::string class_name_ = "ShapeUnion";
  std::shared_ptr<Shape> shape1_, shape2_;
};

inline std::shared_ptr<ShapeUnion> MakeShapeUnion(
    std::shared_ptr<Shape> shape1,
    std::shared_ptr<Shape> shape2) {
  return std::make_shared<ShapeUnion>(shape1, shape2);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_SHAPE_UNION_H_
