
#include "system/include/bond_visitor.h"

namespace feasst {

class MapBondVisitor {
 public:
  MapBondVisitor() {
    auto obj = MakeBondVisitor();
    obj->deserialize_map()["BondVisitor"] = obj;
  }
};

static MapBondVisitor mapper_ = MapBondVisitor();

std::map<std::string, std::shared_ptr<BondVisitor> >& BondVisitor::deserialize_map() {
  static std::map<std::string, std::shared_ptr<BondVisitor> >* ans =
     new std::map<std::string, std::shared_ptr<BondVisitor> >();
  return *ans;
}

void BondVisitor::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bond_visitor_(ostr);
}

std::shared_ptr<BondVisitor> BondVisitor::create(std::istream& istr) const {
  return std::make_shared<BondVisitor>(istr);
}

std::shared_ptr<BondVisitor> BondVisitor::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void BondVisitor::serialize_bond_visitor_(std::ostream& ostr) const {
  feasst_serialize_version(303, ostr);
  feasst_serialize(energy_, ostr);
}

BondVisitor::BondVisitor(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(303 == version, "mismatch version: " << version);
  feasst_deserialize(&energy_, istr);
}

}  // namespace feasst
