#include "configuration/include/site.h"
#include "utils/include/debug.h"

namespace feasst {

void Site::add_property(const std::string name, const double value) {
  PropertiedEntity::add_property(name, value);
  if (name == "director") {
    is_director_ = true;
  }
}

void Site::serialize(std::ostream& ostr) const {
  PropertiedEntity::serialize(ostr);
  TypedEntity::serialize(ostr);
  SpatialEntity::serialize(ostr);
  feasst_serialize_version(1, ostr);
  feasst_serialize(is_director_, ostr);
  feasst_serialize(is_physical_, ostr);
}

Site::Site(std::istream& istr)
  : PropertiedEntity(istr),
    TypedEntity(istr),
    SpatialEntity(istr) {
  feasst_deserialize_version(istr);
  feasst_deserialize(&is_director_, istr);
  feasst_deserialize(&is_physical_, istr);
}

}  // namespace feasst
