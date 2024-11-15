#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "charge/include/check_net_charge.h"

namespace feasst {

CheckNetCharge::CheckNetCharge(argtype args) : AnalyzeUpdateOnly(&args) {
  minimum_ = dble("minimum", &args, 0.);
  maximum_ = dble("maximum", &args, 0.);
  feasst_check_all_used(args);
}

void CheckNetCharge::update(const MonteCarlo& mc) {
  const double net_charge = ewald_.net_charge(configuration(mc.system()));
  ASSERT(net_charge > minimum_ - 100.*NEAR_ZERO &&
         net_charge < maximum_ + 100.*NEAR_ZERO,
    "The net charge: " << net_charge << " should be > " << minimum_ <<
    " and < " << maximum_);
}

FEASST_MAPPER(CheckNetCharge,);

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
