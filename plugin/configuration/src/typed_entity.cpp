
#include "configuration/include/typed_entity.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"

namespace feasst {

void TypedEntity::serialize(std::ostream& ostr) const {
  feasst_serialize_version(8896, ostr);
  feasst_serialize(type_, ostr);
}

TypedEntity::TypedEntity(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8896, "version mismatch: " << version);
  feasst_deserialize(&type_, istr);
}

}  // namespace feasst
