#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "chain/include/select_reptate.h"
#include "chain/include/perturb_reptate.h"
#include "chain/include/trial_reptate.h"

namespace feasst {

FEASST_MAPPER(TrialReptate, argtype({{"max_length", "1"}}));

TrialReptate::TrialReptate(argtype * args) :
  TrialMove(std::make_shared<SelectReptate>(args),
            std::make_shared<PerturbReptate>(args),
            args) {
  class_name_ = "TrialReptate";
  set_description("TrialReptate");
}
TrialReptate::TrialReptate(argtype args) : TrialReptate(&args) {
  feasst_check_all_used(args);
}

TrialReptate::TrialReptate(std::istream& istr) : TrialMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2456, "mismatch version: " << version);
}

void TrialReptate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(2456, ostr);
}

}  // namespace feasst
