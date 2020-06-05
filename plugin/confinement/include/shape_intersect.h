
#ifndef FEASST_CONFINEMENT_SHAPE_INTERSECT_H_
#define FEASST_CONFINEMENT_SHAPE_INTERSECT_H_

#include <memory>
#include <sstream>
#include "confinement/include/shape.h"

namespace feasst {

/**
  Represents the intersection of two shapes.
  An intersection of two shapes is the region which is inside both shapes.
  This is implemented by returning the largest value of the nearest distance.
 */
class ShapeIntersect : public Shape {
 public:
  // This constructor only to be used for serialization.
  ShapeIntersect() {}

  ShapeIntersect(std::shared_ptr<Shape> shape1,
                 std::shared_ptr<Shape> shape2);

  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<ShapeIntersect>(istr); }
  explicit ShapeIntersect(std::istream& istr);
  virtual ~ShapeIntersect() {}

 private:
  const std::string class_name_ = "ShapeIntersect";
  std::shared_ptr<Shape> shape1_, shape2_;
};

inline std::shared_ptr<ShapeIntersect> MakeShapeIntersect() {
  return std::make_shared<ShapeIntersect>();
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_SHAPE_INTERSECT_H_
