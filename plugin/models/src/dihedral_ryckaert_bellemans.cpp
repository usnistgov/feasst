#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "configuration/include/bond.h"
#include "models/include/dihedral_ryckaert_bellemans.h"

namespace feasst {

class MapDihedralRyckaertBellemans {
 public:
  MapDihedralRyckaertBellemans() {
    auto obj = MakeDihedralRyckaertBellemans();
    obj->deserialize_map()["DihedralRyckaertBellemans"] = obj;
  }
};

static MapDihedralRyckaertBellemans mapper_ = MapDihedralRyckaertBellemans();

std::shared_ptr<BondFourBody> DihedralRyckaertBellemans::create(std::istream& istr) const {
  return std::make_shared<DihedralRyckaertBellemans>(istr);
}

DihedralRyckaertBellemans::DihedralRyckaertBellemans(std::istream& istr) : BondFourBody(istr) {
  // ASSERT(class_name_ == "DihedralRyckaertBellemans", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7048 == version, "mismatch version: " << version);
}

void DihedralRyckaertBellemans::serialize_dihedral_ryckaert_bellemans_(std::ostream& ostr) const {
  serialize_bond_four_body_(ostr);
  feasst_serialize_version(7048, ostr);
}

void DihedralRyckaertBellemans::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_dihedral_ryckaert_bellemans_(ostr);
}

double DihedralRyckaertBellemans::energy(const double radians, const Bond& dihedral) const {
  const double c0 = dihedral.property("c0");
  const double c1 = dihedral.property("c1");
  const double c2 = dihedral.property("c2");
  const double c3 = dihedral.property("c3");
  const double cosphi = std::cos(radians);
  const double en = c0 + cosphi*(c1 + cosphi*(c2 + cosphi*c3));
  ASSERT(!std::isnan(en) && !std::isinf(en), "en: " << en << " radians: "
    << radians << " c0 " << c0 << " c1 " << c1 << " c2 " << c2 << " c3 " << c3);
  return en;
}

}  // namespace feasst
