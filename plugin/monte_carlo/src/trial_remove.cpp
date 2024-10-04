#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/trial_remove.h"

namespace feasst {

FEASST_MAPPER(TrialRemove,);

TrialRemove::TrialRemove(argtype * args) : Trial(args) {
  class_name_ = "TrialRemove";
  set_description("TrialRemove");
  // optimization: do not load coordinates if num_steps == 1, by default
  { argtype tmpargs = *args;
    if (integer("num_steps", &tmpargs, 1) == 1) {
      if (!used("load_coordinates", tmpargs)) {
        args->insert({"load_coordinates", "false"});
      }
    }
  }
  add_stage(
    std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbRemove>(),
    args
  );
  set(std::make_shared<TrialComputeRemove>(args));
}
TrialRemove::TrialRemove(argtype args) : TrialRemove(&args) {
  feasst_check_all_used(args);
}

TrialRemove::TrialRemove(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2304, "mismatch version: " << version);
}

void TrialRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(2304, ostr);
}

}  // namespace feasst
