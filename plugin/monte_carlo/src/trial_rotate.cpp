#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "monte_carlo/include/trial_rotate.h"

namespace feasst {

FEASST_MAPPER(TrialRotate,);

TrialRotate::TrialRotate(argtype * args) :
  TrialMove(std::make_shared<TrialSelectParticle>(args),
            std::make_shared<PerturbRotate>(args),
            args) {
  class_name_ = "TrialRotate";
  set_description("TrialRotate");
}
TrialRotate::TrialRotate(argtype args) : TrialRotate(&args) {
  feasst_check_all_used(args);
}

TrialRotate::TrialRotate(std::istream& istr) : TrialMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1867, "mismatch version: " << version);
}

void TrialRotate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(1867, ostr);
}

}  // namespace feasst
