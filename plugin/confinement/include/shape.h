
#ifndef FEASST_CONFINEMENT_SHAPE_H_
#define FEASST_CONFINEMENT_SHAPE_H_

#include <memory>
#include <map>
#include <sstream>
#include "core/include/position.h"
#include "core/include/utils_io.h"

namespace feasst {

/**
  Shapes may be defined by either a simple mathematical
  formula or an interpolated data table.
 */
class Shape {
 public:
  /// Return the distance from the point to the nearest point on the surface.
  /// The distance is negative if the point is inside of the shape.
  virtual double nearest_distance(const Position& point) const = 0;

  /// Return true if the point is inside of the shape.
  bool is_inside(const Position& point) const;

  /// Return true if the sphere of given center point and diameter is entirely
  /// inside of the shape.
  bool is_inside(const Position& point, const double diameter) const;

  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Shape> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Shape> >& deserialize_map();
  std::shared_ptr<Shape> deserialize(std::istream& istr);
  virtual ~Shape() {}
};

// An object which contains a shape.
class ShapedEntity {
 public:
  ShapedEntity() {}
  ShapedEntity(std::shared_ptr<Shape> shape) {
    shape_ = shape; }

  /// Set the shape.
  void set_shape(const std::shared_ptr<Shape> shape) { shape_ = shape; }

  /// Return the shape.
  const std::shared_ptr<Shape> shape() const { return shape_; }

  void serialize(std::ostream& ostr) const {
    feasst_serialize_version(1, ostr);
    feasst_serialize_fstdr(shape_, ostr);
  }

  ShapedEntity(std::istream& istr) {
    feasst_deserialize_version(istr);
    // feasst_deserialize_fstdr(shape_, istr);
    { // HWH for unknown reasons the above template function does not work
      int existing;
      istr >> existing;
      if (existing != 0) {
        shape_ = shape_->deserialize(istr);
      }
    }
  }

 private:
  std::shared_ptr<Shape> shape_;
};

/**
  Represents the intersection of two shapes.
  This is done by returning the largest value of the nearest distance.
 */
class ShapeIntersect : public Shape {
 public:
  // This constructor only to be used for serialization.
  ShapeIntersect() {}

  ShapeIntersect(std::shared_ptr<Shape> shape1,
                 std::shared_ptr<Shape> shape2);

  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(822, ostr);
    feasst_serialize_fstdr(shape1_, ostr);
    feasst_serialize_fstdr(shape2_, ostr);
  }

  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<ShapeIntersect>(istr); }

  ShapeIntersect(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(822 == version, version);

    // HWH for unknown reasons, below template isn't working in this case.
    // feasst_deserialize_fstdr(shape1, istr);
    // feasst_deserialize_fstdr(shape2, istr);
    int existing;
    istr >> existing;
    if (existing != 0) {
      shape1_ = shape1_->deserialize(istr);
    }
    istr >> existing;
    if (existing != 0) {
      shape2_ = shape2_->deserialize(istr);
    }
  }

  virtual ~ShapeIntersect() {}

 private:
  const std::string class_name_ = "ShapeIntersect";
  std::shared_ptr<Shape> shape1_, shape2_;
};

inline std::shared_ptr<ShapeIntersect> MakeShapeIntersect() {
  return std::make_shared<ShapeIntersect>();
}

// HWH encompase Intersect and Union with a single two-shape base class.
/**
  Represents the union of two shapes.
  This is done by returning the smallest value of the nearest distance.
 */
class ShapeUnion : public Shape {
 public:
  // This constructor only to be used for serialization.
  ShapeUnion() {}

  ShapeUnion(std::shared_ptr<Shape> shape1,
             std::shared_ptr<Shape> shape2);

  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(172, ostr);
    feasst_serialize_fstdr(shape1_, ostr);
    feasst_serialize_fstdr(shape2_, ostr);
  }

  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<ShapeUnion>(istr); }

  ShapeUnion(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(172 == version, version);

    // HWH for unknown reasons, below template isn't working in this case.
    // feasst_deserialize_fstdr(shape1, istr);
    // feasst_deserialize_fstdr(shape2, istr);
    int existing;
    istr >> existing;
    if (existing != 0) {
      shape1_ = shape1_->deserialize(istr);
    }
    istr >> existing;
    if (existing != 0) {
      shape2_ = shape2_->deserialize(istr);
    }
  }

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

#endif  // FEASST_CONFINEMENT_SHAPE_H_
