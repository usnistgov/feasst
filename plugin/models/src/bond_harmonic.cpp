#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "configuration/include/bond.h"
#include "models/include/bond_harmonic.h"

namespace feasst {

class MapBondHarmonic {
 public:
  MapBondHarmonic() {
    auto obj = MakeBondHarmonic();
    obj->deserialize_map()["BondHarmonic"] = obj;
  }
};

static MapBondHarmonic mapper_ = MapBondHarmonic();

std::shared_ptr<BondTwoBody> BondHarmonic::create(std::istream& istr) const {
  return std::make_shared<BondHarmonic>(istr);
}

BondHarmonic::BondHarmonic(std::istream& istr) : BondTwoBody(istr) {
  // ASSERT(class_name_ == "BondHarmonic", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(3864 == version, "mismatch version: " << version);
}

void BondHarmonic::serialize_bond_harmonic_(std::ostream& ostr) const {
  serialize_bond_two_body_(ostr);
  feasst_serialize_version(3864, ostr);
}

void BondHarmonic::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bond_harmonic_(ostr);
}

double BondHarmonic::energy(const double distance, const Bond& bond) const {
  const double k = bond.property("k_energy_per_length_sq");
  const double l0 = bond.property("equilibrium_length");
  const double dl = distance - l0;
  const double en = k*dl*dl;
  DEBUG("bond harmonic " << en);
  return en;
}

double BondHarmonic::random_distance(const Bond& bond, const double beta,
    const int dimen, Random * random) const {
  const double equilibrium_length = bond.property("equilibrium_length");
  const double beta_k = beta*bond.property("k_energy_per_length_sq");
  const double sigma = std::sqrt(1./2./beta_k);
  const double max_length_dm1 = std::pow(equilibrium_length + 3.*sigma, dimen - 1);
  int attempt = 0;
  while (attempt < 1e6) {
    if (dimen == 3) {
      const double length = random->normal(equilibrium_length, sigma);
      if (random->uniform() < length*length/max_length_dm1) return length;
    } else if (dimen == 2) {
      const double length = random->uniform_real(equilibrium_length - 3.*sigma,
                                                 equilibrium_length + 3.*sigma);
      const double dl = length - equilibrium_length;
      if (random->uniform() < (length/max_length_dm1)*std::exp(-beta_k*dl*dl)) return length;
    }
    ++attempt;
  }
  FATAL("maximum attempts reached");
}

}  // namespace feasst
