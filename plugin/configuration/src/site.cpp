#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "configuration/include/site.h"

namespace feasst {

Site::Site() : PropertiedEntity(), TypedEntity() {
  set_physical();
  set_anisotropic();
  is_anisotropic_ = false;
}

void Site::add_property(const std::string name, const double value) {
  PropertiedEntity::add_property(name, value);
}

int Site::cell(const int index) const {
  ASSERT(index < num_cells(),
    "index: " << index << " cellsize: " << cells_.size());
  return cells_[index];
}

void Site::set_euler(const Euler& euler) {
  euler_ = euler;
  is_anisotropic_ = true;
}

void Site::set_anisotropic(const bool aniso) {
  is_anisotropic_ = aniso;
}

bool Site::is_anisotropic() const {
  return is_anisotropic_;
}

void Site::serialize(std::ostream& ostr) const {
  PropertiedEntity::serialize(ostr);
  TypedEntity::serialize(ostr);
  feasst_serialize_version(481, ostr);
  feasst_serialize_fstobj(position_, ostr);
  feasst_serialize_fstobj(euler_, ostr);
  feasst_serialize(is_physical_, ostr);
  feasst_serialize(is_anisotropic_, ostr);
  feasst_serialize(cells_, ostr);
}

Site::Site(std::istream& istr)
  : PropertiedEntity(istr),
    TypedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 480 || version == 481, "unrecognized version: " << version);
  feasst_deserialize_fstobj(&position_, istr);
  feasst_deserialize_fstobj(&euler_, istr);
  if (version == 480) {
    bool tmp;
    feasst_deserialize(&tmp, istr);
  }
  feasst_deserialize(&is_physical_, istr);
  feasst_deserialize(&is_anisotropic_, istr);
  feasst_deserialize(&cells_, istr);
}

}  // namespace feasst
