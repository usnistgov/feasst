#include <cmath>
#include "utils/include/serialize.h"
#include "chain/include/analyze_rigid_bonds.h"
#include "math/include/constants.h"

namespace feasst {

class MapAnalyzeRigidBonds {
 public:
  MapAnalyzeRigidBonds() {
    AnalyzeRigidBonds().deserialize_map()["AnalyzeRigidBonds"] =
      MakeAnalyzeRigidBonds();
  }
};

static MapAnalyzeRigidBonds mapper_ = MapAnalyzeRigidBonds();

AnalyzeRigidBonds::AnalyzeRigidBonds(const argtype &args) : AnalyzeUpdateOnly(args) {}

void AnalyzeRigidBonds::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(549, ostr);
  feasst_serialize_fstobj(visitor_, ostr);
  feasst_serialize_fstobj(bond_, ostr);
  feasst_serialize_fstobj(angle_, ostr);
}

AnalyzeRigidBonds::AnalyzeRigidBonds(std::istream& istr)
  : AnalyzeUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 549, "version mismatch: " << version);
  feasst_deserialize_fstobj(&visitor_, istr);
  feasst_deserialize_fstobj(&bond_, istr);
  feasst_deserialize_fstobj(&angle_, istr);
}

void AnalyzeRigidBonds::update(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  visitor_.compute(bond_, system.configuration());
  ASSERT(std::abs(visitor_.energy()) < NEAR_ZERO, "bond check failure");
  visitor_.compute(angle_, system.configuration());
  ASSERT(std::abs(visitor_.energy()) < NEAR_ZERO, "angle check failure");
}

}  // namespace feasst
