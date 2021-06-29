#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

TrialMove::TrialMove(std::shared_ptr<TrialSelect> select,
  std::shared_ptr<PerturbMove> perturb,
  argtype * args) : Trial(args) {
  add_stage(select, perturb, args);
  set(std::make_shared<TrialComputeMove>(args));
}

TrialMove::TrialMove(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3294, "mismatch version: " << version);
}

void TrialMove::serialize_trial_move_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(3294, ostr);
}

}  // namespace feasst
