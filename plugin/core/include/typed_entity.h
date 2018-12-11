
#ifndef FEASST_CORE_TYPED_ENTITY_H_
#define FEASST_CORE_TYPED_ENTITY_H_

#include <vector>
#include <string>

namespace feasst {

/**
  Manage entity types.
 */
class TypedEntity {
 public:
  /// Obtain the type.
  int type() const { return type_; }

  /// Set the type.
  void set_type(const int type) { type_ = type; }

 private:
  int type_ = 0;
};

}  // namespace feasst

#endif  // FEASST_CORE_TYPED_ENTITY_H_
