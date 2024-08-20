
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "configuration/include/bond.h"

namespace feasst {

Bond::Bond(std::istream& istr) : PropertiedEntity(istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 743, "version mismatch: " << version);
  feasst_deserialize(&type_, istr);
  feasst_deserialize(&site_indicies_, istr);
  feasst_deserialize(&model_, istr);
}

void Bond::serialize(std::ostream& ostr) const {
  PropertiedEntity::serialize(ostr);
  ostr << class_name_ << " ";
  feasst_serialize_version(743, ostr);
  feasst_serialize(type_, ostr);
  feasst_serialize(site_indicies_, ostr);
  feasst_serialize(model_, ostr);
}

}  // namespace feasst
