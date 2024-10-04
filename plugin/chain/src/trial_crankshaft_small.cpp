#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/position.h"
#include "chain/include/select_crankshaft_small.h"
#include "chain/include/perturb_crankshaft_small.h"
#include "chain/include/trial_crankshaft_small.h"

namespace feasst {

FEASST_MAPPER(TrialCrankshaftSmall,
  argtype({{"site", "0"}, {"anchor_site0", "1"}, {"anchor_site1", "2"}}));

TrialCrankshaftSmall::TrialCrankshaftSmall(argtype * args) :
  TrialMove(std::make_shared<SelectCrankshaftSmall>(args),
            std::make_shared<PerturbCrankshaftSmall>(args),
            args) {
  class_name_ = "TrialCrankshaftSmall";
  set_description("TrialCrankshaftSmall");
}
TrialCrankshaftSmall::TrialCrankshaftSmall(argtype args) : TrialCrankshaftSmall(&args) {
  feasst_check_all_used(args);
}

TrialCrankshaftSmall::TrialCrankshaftSmall(std::istream& istr) : TrialMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1826, "mismatch version: " << version);
}

void TrialCrankshaftSmall::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(1826, ostr);
}

}  // namespace feasst
