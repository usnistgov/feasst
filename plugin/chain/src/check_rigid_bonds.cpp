#include <cmath>
#include "utils/include/serialize.h"
#include "chain/include/check_rigid_bonds.h"
#include "math/include/constants.h"

namespace feasst {

class MapCheckRigidBonds {
 public:
  MapCheckRigidBonds() {
    CheckRigidBonds().deserialize_map()["CheckRigidBonds"] =
      MakeCheckRigidBonds();
  }
};

static MapCheckRigidBonds mapper_ = MapCheckRigidBonds();

CheckRigidBonds::CheckRigidBonds(const argtype &args)
  : AnalyzeUpdateOnly(args) {
  visitor_ = MakeBondVisitor({{"verbose", "true"}});
}

void CheckRigidBonds::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(549, ostr);
  feasst_serialize(visitor_, ostr);
  feasst_serialize_fstobj(bond_, ostr);
  feasst_serialize_fstobj(angle_, ostr);
}

CheckRigidBonds::CheckRigidBonds(std::istream& istr)
  : AnalyzeUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 549, "version mismatch: " << version);
  // HWH for unknown reasons, this function template does not work.
  /// Deserialize feasst object stored as shared pointer
  //feasst_deserialize(visitor_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      visitor_ = std::make_shared<BondVisitor>(istr);
    }
  }

  feasst_deserialize_fstobj(&bond_, istr);
  feasst_deserialize_fstobj(&angle_, istr);
}

void CheckRigidBonds::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  visitor_->compute(bond_, system.configuration());
  ASSERT(std::abs(visitor_->energy()) < NEAR_ZERO, "bond check failure");
  visitor_->compute(angle_, system.configuration());
  ASSERT(std::abs(visitor_->energy()) < NEAR_ZERO, "angle check failure");
}

}  // namespace feasst
