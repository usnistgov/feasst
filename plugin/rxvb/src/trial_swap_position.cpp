#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "rxvb/include/select_two_particles.h"
#include "rxvb/include/perturb_swap_position.h"
#include "rxvb/include/trial_swap_position.h"

namespace feasst {

FEASST_MAPPER(TrialSwapPosition,);

TrialSwapPosition::TrialSwapPosition(argtype * args) :
  TrialMove(std::make_shared<SelectTwoParticles>(args),
            std::make_shared<PerturbSwapPosition>(args),
            args) {
  class_name_ = "TrialSwapPosition";
  set_description("TrialSwapPosition");
}
TrialSwapPosition::TrialSwapPosition(argtype args) : TrialSwapPosition(&args) {
  feasst_check_all_used(args);
}

TrialSwapPosition::TrialSwapPosition(std::istream& istr) : TrialMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6705, "mismatch version: " << version);
}

void TrialSwapPosition::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(6705, ostr);
}

}  // namespace feasst
