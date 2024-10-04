#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "system/include/system.h"
#include "gibbs/include/check_constant_volume.h"

namespace feasst {

FEASST_MAPPER(CheckConstantVolume,);

CheckConstantVolume::CheckConstantVolume(argtype * args) : ModifyUpdateOnly(args) {
  tolerance_ = dble("tolerance", args, 1e-4);
}
CheckConstantVolume::CheckConstantVolume(argtype args) : CheckConstantVolume(&args) {
  feasst_check_all_used(args);
}

void CheckConstantVolume::update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) {
  const double total_volume = system->total_volume();
  if (last_total_volume_ < 0) {
    last_total_volume_ = total_volume;
    return;
  }
  const double diff = total_volume - last_total_volume_;
  ASSERT(std::abs(diff) < tolerance_,
    "The difference in total volume, " << diff << " is greater than the " <<
    "tolerance: " << tolerance_ << " . The last recorded volume was " <<
    last_total_volume_ << " while the current is " << total_volume);
  last_total_volume_ = total_volume;
}

void CheckConstantVolume::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(7024, ostr);
  feasst_serialize(tolerance_, ostr);
}

CheckConstantVolume::CheckConstantVolume(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7024, "version mismatch: " << version);
  feasst_deserialize(&tolerance_, istr);
}

}  // namespace feasst
