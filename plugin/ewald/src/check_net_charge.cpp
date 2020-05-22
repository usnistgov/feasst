#include "ewald/include/check_net_charge.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"

namespace feasst {

CheckNetCharge::CheckNetCharge(const argtype &args) : AnalyzeUpdateOnly(args) {
  minimum_ = args_.key("minimum").dflt("0.").dble();
  maximum_ = args_.key("maximum").dflt("0.").dble();
}

void CheckNetCharge::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  const double net_charge = ewald_.net_charge(system.configuration());
  ASSERT(net_charge > minimum_ - 100.*NEAR_ZERO &&
         net_charge < maximum_ + 100.*NEAR_ZERO,
    "The net charge: " << net_charge << " should be > " << minimum_ <<
    " and < " << maximum_);
}

class MapCheckNetCharge {
 public:
  MapCheckNetCharge() {
    CheckNetCharge().deserialize_map()["CheckNetCharge"] =
      MakeCheckNetCharge();
  }
};

static MapCheckNetCharge mapper_ = MapCheckNetCharge();

void CheckNetCharge::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(943, ostr);
  feasst_serialize(minimum_, ostr);
  feasst_serialize(maximum_, ostr);
}

CheckNetCharge::CheckNetCharge(std::istream& istr) : AnalyzeUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 943, "version mismatch: " << version);
  feasst_deserialize(&minimum_, istr);
  feasst_deserialize(&maximum_, istr);
}

}  // namespace feasst
