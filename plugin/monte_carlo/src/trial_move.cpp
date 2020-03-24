#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

TrialMove::TrialMove(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialMove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(691 == version, "mismatch version: " << version);
}

void TrialMove::serialize_trial_move_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(691, ostr);
}

}  // namespace feasst
