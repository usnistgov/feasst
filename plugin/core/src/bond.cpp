
#include "core/include/bond.h"
#include "core/include/debug.h"
#include "core/include/utils_io.h"

namespace feasst {

Bond::Bond(std::istream& istr)
  : PropertiedEntity(istr),
    TypedEntity(istr) {
  int version;
  istr >> version;
  istr >> name_;
  feasst_deserialize(&site_indicies_, istr);
}

void Bond::serialize(std::ostream& ostr) const {
  PropertiedEntity::serialize(ostr);
  TypedEntity::serialize(ostr);
  ostr << MAX_PRECISION;
  ostr << "1 "; // version
  ostr << name_ << " ";
  feasst_serialize(site_indicies_, ostr);
}

}  // namespace feasst
