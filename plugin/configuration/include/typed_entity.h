
#ifndef FEASST_CONFIGURATION_TYPED_ENTITY_H_
#define FEASST_CONFIGURATION_TYPED_ENTITY_H_

#include <vector>
#include <string>

namespace feasst {

/**
  Manage entity types.
 */
class TypedEntity {
 public:
  TypedEntity() {}

  /// Obtain the type.
  int type() const { return type_; }

  /// Set the type.
  void set_type(const int type) { type_ = type; }

  void serialize(std::ostream& ostr) const;
  TypedEntity(std::istream& istr);

 private:
  int type_ = 0;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_TYPED_ENTITY_H_
