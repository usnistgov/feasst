#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "models/include/fene.h"

namespace feasst {

class MapFENE {
 public:
  MapFENE() {
    auto obj = MakeFENE();
    obj->deserialize_map()["FENE"] = obj;
  }
};

static MapFENE mapper_ = MapFENE();

std::shared_ptr<BondTwoBody> FENE::create(std::istream& istr) const {
  return std::make_shared<FENE>(istr);
}

FENE::FENE(std::istream& istr) : BondTwoBody(istr) {
  // ASSERT(class_name_ == "FENE", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(2945 == version, "mismatch version: " << version);
}

void FENE::serialize_fene_(std::ostream& ostr) const {
  serialize_bond_two_body_(ostr);
  feasst_serialize_version(2945, ostr);
}

void FENE::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_fene_(ostr);
}

double FENE::energy(
    const double distance,
    const Bond& bond) const {
  const double R0 = bond.property("R0");
  if (distance >= R0) {
    return NEAR_INFINITY;
  }
  const double k = bond.property("k_energy_per_length_sq");
  return energy_fene(distance, k, R0);
}

double FENE::energy_fene(const double distance, const double k,
    const double R0) const {
  TRACE("distance " << distance);
  TRACE("k " << k);
  TRACE("R0 " << R0);
  const double en = -0.5*k*R0*R0*std::log(1.-distance*distance/R0/R0);
  TRACE("en " << en);
  return en;
}

double FENE::random_distance(const Bond& bond, const double beta,
    const int dimen, Random * random) const {
  const double R0 = bond.property("R0");
  const double k = bond.property("k_energy_per_length_sq");
  double jacobian = 0;
  const double max_length_dm1 = std::pow(R0, dimen - 1);
  int attempt = 0;
  while (attempt < 1e6) {
    const double length = random->uniform_real(0, R0);
    if (dimen == 3) {
      jacobian = length*length;
    } else if (dimen == 2) {
      jacobian = length;
    }
    if (random->uniform() <
        jacobian/max_length_dm1*std::exp(-beta*energy_fene(length, k, R0))) {
      return length;
    }
    ++attempt;
  }
  FATAL("maximum attempts reached");
}

}  // namespace feasst
