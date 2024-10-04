#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/trial_translate.h"

namespace feasst {

FEASST_MAPPER(TrialTranslate,);

TrialTranslate::TrialTranslate(argtype * args) :
  TrialMove(std::make_shared<TrialSelectParticle>(args),
            std::make_shared<PerturbTranslate>(args),
            args) {
  class_name_ = "TrialTranslate";
  set_description("TrialTranslate");
}
TrialTranslate::TrialTranslate(argtype args) : TrialTranslate(&args) {
  feasst_check_all_used(args);
}

TrialTranslate::TrialTranslate(std::istream& istr) : TrialMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3056, "mismatch version: " << version);
}

void TrialTranslate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(3056, ostr);
}

}  // namespace feasst
