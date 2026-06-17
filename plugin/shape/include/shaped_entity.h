
#ifndef FEASST_SHAPE_SHAPED_ENTITY_H_
#define FEASST_SHAPE_SHAPED_ENTITY_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace feasst {

class Shape;

/**
  An object which contains a shape.
 */
class ShapedEntity {
 public:
  ShapedEntity() {}
  explicit ShapedEntity(std::shared_ptr<Shape> shape);

  void set_shape(std::shared_ptr<Shape> shape);

  /// Return the shape.
  const std::shared_ptr<Shape> shape() const;

  void serialize(std::ostream& ostr) const;
  explicit ShapedEntity(std::istream& istr);

 private:
  std::shared_ptr<Shape> shape_;
};

}  // namespace feasst

#endif  // FEASST_SHAPE_SHAPED_ENTITY_H_
