
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "system/include/bond_three_body.h"

namespace feasst {

std::map<std::string, std::shared_ptr<BondThreeBody> >& BondThreeBody::deserialize_map() {
  static std::map<std::string, std::shared_ptr<BondThreeBody> >* ans =
     new std::map<std::string, std::shared_ptr<BondThreeBody> >();
  return *ans;
}

void BondThreeBody::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<BondThreeBody> BondThreeBody::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<BondThreeBody> BondThreeBody::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void BondThreeBody::serialize_bond_three_body_(std::ostream& ostr) const {
  feasst_serialize_version(943, ostr);
}

BondThreeBody::BondThreeBody(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(943 == version, "mismatch version: " << version);
}

}  // namespace feasst
