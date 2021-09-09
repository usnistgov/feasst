#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "models/include/dihedral_harmonic.h"

namespace feasst {

class MapDihedralHarmonic {
 public:
  MapDihedralHarmonic() {
    auto obj = MakeDihedralHarmonic();
    obj->deserialize_map()["DihedralHarmonic"] = obj;
  }
};

static MapDihedralHarmonic mapper_ = MapDihedralHarmonic();

std::shared_ptr<BondFourBody> DihedralHarmonic::create(std::istream& istr) const {
  return std::make_shared<DihedralHarmonic>(istr);
}

DihedralHarmonic::DihedralHarmonic(std::istream& istr) : BondFourBody(istr) {
  // ASSERT(class_name_ == "DihedralHarmonic", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(8385 == version, "mismatch version: " << version);
}

void DihedralHarmonic::serialize_dihedral_harmonic_(std::ostream& ostr) const {
  serialize_bond_four_body_(ostr);
  feasst_serialize_version(8385, ostr);
}

void DihedralHarmonic::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_dihedral_harmonic_(ostr);
}

double DihedralHarmonic::energy(const double radians, const Bond& dihedral) const {
  DEBUG("radians " << radians);
  const double equil_radians =
    degrees_to_radians(dihedral.property("equilibrium_degrees"));
  const double k = dihedral.property("k_energy_per_radian_sq");
  double delta_rad = radians - equil_radians;
  DEBUG("delta_rad " << delta_rad);
  return k*delta_rad*delta_rad;
}

}  // namespace feasst
