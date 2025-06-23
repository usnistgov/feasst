#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add_remove.h"
#include "monte_carlo/include/trial_compute_add_remove.h"
#include "monte_carlo/include/trial_add_remove.h"

namespace feasst {

FEASST_MAPPER(TrialAddRemove,);

TrialAddRemove::TrialAddRemove(argtype * args) : Trial(args) {
  class_name_ = "TrialAddRemove";
  set_description("TrialAddRemove");
  auto perturb = std::make_shared<PerturbAddRemove>(args);
  args->insert({"half_ghost", "true"});
  add_stage(std::make_shared<TrialSelectParticle>(args), perturb, args);
  set(std::make_shared<TrialComputeAddRemove>(args));
}
TrialAddRemove::TrialAddRemove(argtype args) : TrialAddRemove(&args) {
  feasst_check_all_used(args);
}
TrialAddRemove::~TrialAddRemove() {}

TrialAddRemove::TrialAddRemove(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3986, "mismatch version: " << version);
}

void TrialAddRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(3986, ostr);
}

}  // namespace feasst
