
#include "core/include/typed_entity.h"
#include "core/include/debug.h"

namespace feasst {

void TypedEntity::serialize(std::ostream& ostr) const {
  ostr << MAX_PRECISION;
  ostr << "1 "; // version
  ostr << type_ << " ";
}

TypedEntity::TypedEntity(std::istream& istr) {
  int version;
  istr >> version;
  istr >> type_;
}

}  // namespace feasst
