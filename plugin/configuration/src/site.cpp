#include "configuration/include/site.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"

namespace feasst {

Site::Site() : PropertiedEntity(), TypedEntity() {
  set_physical();
}

void Site::add_property(const std::string name, const double value) {
  PropertiedEntity::add_property(name, value);
  if (name == "director") {
    is_director_ = true;
  }
}

int Site::cell(const int index) const {
  ASSERT(index < num_cells(),
    "index: " << index << " cellsize: " << cells_.size());
  return cells_[index];
}

void Site::serialize(std::ostream& ostr) const {
  PropertiedEntity::serialize(ostr);
  TypedEntity::serialize(ostr);
  feasst_serialize_version(480, ostr);
  feasst_serialize_fstobj(position_, ostr);
  feasst_serialize(is_director_, ostr);
  feasst_serialize(is_physical_, ostr);
  feasst_serialize(cells_, ostr);
}

Site::Site(std::istream& istr)
  : PropertiedEntity(istr),
    TypedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 480, "unrecognized version: " << version);
  feasst_deserialize_fstobj(&position_, istr);
  feasst_deserialize(&is_director_, istr);
  feasst_deserialize(&is_physical_, istr);
  feasst_deserialize(&cells_, istr);
}

}  // namespace feasst
