
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "system/include/bond_four_body.h"

namespace feasst {

std::map<std::string, std::shared_ptr<BondFourBody> >& BondFourBody::deserialize_map() {
  static std::map<std::string, std::shared_ptr<BondFourBody> >* ans =
     new std::map<std::string, std::shared_ptr<BondFourBody> >();
  return *ans;
}

void BondFourBody::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<BondFourBody> BondFourBody::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<BondFourBody> BondFourBody::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void BondFourBody::serialize_bond_four_body_(std::ostream& ostr) const {
  feasst_serialize_version(7509, ostr);
}

BondFourBody::BondFourBody(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(7509 == version, "mismatch version: " << version);
}

double BondFourBody::energy(const Position& ri, const Position& rj,
    const Position& rk, const Position& rl, const Dihedral& dihedral) const {
  return energy(ri.torsion_angle_radians(rj, rk, rl), dihedral);
}

double BondFourBody::random_dihedral_radians(const Dihedral& dihedral,
    const double beta, const int dimension, Random * random) const {
  ASSERT(dimension == 3, "dihedrals only implemented in 3D");
  int attempt = 0;
  while (attempt < 1e6) {
    const double radians = 2*PI*random->uniform();
    const double en = energy(radians, dihedral);
    if (random->uniform() < std::exp(-beta*en)) {
      return radians;
    }
    ++attempt;
  }
  FATAL("max attempts reached");
}

}  // namespace feasst
