#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "models/include/dihedral_trappe.h"

namespace feasst {

class MapDihedralTraPPE {
 public:
  MapDihedralTraPPE() {
    auto obj = MakeDihedralTraPPE();
    obj->deserialize_map()["DihedralTraPPE"] = obj;
  }
};

static MapDihedralTraPPE mapper_ = MapDihedralTraPPE();

std::shared_ptr<BondFourBody> DihedralTraPPE::create(std::istream& istr) const {
  return std::make_shared<DihedralTraPPE>(istr);
}

DihedralTraPPE::DihedralTraPPE(std::istream& istr) : BondFourBody(istr) {
  // ASSERT(class_name_ == "DihedralTraPPE", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(846 == version, "mismatch version: " << version);
}

void DihedralTraPPE::serialize_dihedral_trappe_(std::ostream& ostr) const {
  serialize_bond_four_body_(ostr);
  feasst_serialize_version(846, ostr);
}

void DihedralTraPPE::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_dihedral_trappe_(ostr);
}

double DihedralTraPPE::energy(const double radians, const Bond& dihedral) const {
  const double c0 = dihedral.property("c0");
  const double c1 = dihedral.property("c1");
  const double c2 = dihedral.property("c2");
  const double c3 = dihedral.property("c3");
  const double en = c0 + c1*(1. + std::cos(radians))
                       + c2*(1. - std::cos(2.*radians))
                       + c3*(1. + std::cos(3.*radians));
  ASSERT(!std::isnan(en) && !std::isinf(en), "en: " << en << " radians: "
    << radians << " c0 " << c0 << " c1 " << c1 << " c2 " << c2 << " c3 " << c3);
  return en;
}

}  // namespace feasst
