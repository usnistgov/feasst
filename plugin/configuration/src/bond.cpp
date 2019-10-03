
#include "configuration/include/bond.h"
#include "utils/include/debug.h"
#include "utils/include/utils_io.h"

namespace feasst {

Bond::Bond(std::istream& istr)
  : PropertiedEntity(istr),
    TypedEntity(istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 743, "version mismatch: " << version);
  feasst_deserialize(&site_indicies_, istr);
}

void Bond::serialize(std::ostream& ostr) const {
  PropertiedEntity::serialize(ostr);
  TypedEntity::serialize(ostr);
  ostr << class_name_ << " ";
  feasst_serialize_version(743, ostr);
  feasst_serialize(site_indicies_, ostr);
}

}  // namespace feasst
