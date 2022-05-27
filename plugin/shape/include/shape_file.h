
#ifndef FEASST_SHAPE_SHAPE_FILE_H_
#define FEASST_SHAPE_SHAPE_FILE_H_

#include <memory>
#include <sstream>
#include "shape/include/shape.h"

namespace feasst {

/**
  Build a Shape using a file with the following syntax.
  The first line describes a Shape by first listing its class_name,
  then its arguements, space-separated.
  The following lines then begin with either "union" or "intersect" and then
  another shape description as described for the first line.
 */
class ShapeFile : public Shape {
 public:
  /**
    args:
    - file_name: name of the file which describes the shape.
   */
  explicit ShapeFile(argtype args);
  explicit ShapeFile(argtype * args);

  bool is_inside(const Position& point) const override {
    return shape_->is_inside(point); }
  bool is_inside(const Position& point, const double diameter) const override {
    return shape_->is_inside(point, diameter); }
  double nearest_distance(const Position& point) const override {
    return nearest_distance(point); }

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<ShapeFile>(istr); }
  explicit ShapeFile(std::istream& istr);
  virtual ~ShapeFile() {}

 private:
  std::shared_ptr<Shape> shape_;
};

inline std::shared_ptr<ShapeFile> MakeShapeFile(argtype args) {
  return std::make_shared<ShapeFile>(args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_SHAPE_FILE_H_
