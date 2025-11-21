#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "morph/include/perturb_particle_type.h"
#include "morph/include/trial_position_swap.h"

namespace feasst {

FEASST_MAPPER(TrialPositionSwap, argtype({
  {"particle_type", "0"}, {"particle_type_morph", "1"}}));

TrialPositionSwap::TrialPositionSwap(argtype * args) : Trial(args) {
  class_name_ = "TrialPositionSwap";
  set_description("TrialPositionSwap");
  DEBUG("args: " << str(*args));
  const std::string t1 = str("particle_type", args);
  const std::string t2 = str("particle_type_morph", args);
  DEBUG("t1: " << t1 << " t2: " << t2);
  argtype first_args = *args;
  first_args.insert({"particle_type", t1});
  add_stage(
    std::make_shared<TrialSelectParticle>(&first_args),
    std::make_shared<PerturbParticleType>(argtype({{"type", t2}})),
    &first_args);
  feasst_check_all_used(first_args);
  args->insert({"particle_type", t2});
  add_stage(
    std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbParticleType>(argtype({{"type", t1}})),
    args);
  set(std::make_shared<TrialComputeMove>());
}
TrialPositionSwap::TrialPositionSwap(argtype args) : TrialPositionSwap(&args) {
  feasst_check_all_used(args);
}

TrialPositionSwap::TrialPositionSwap(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7934, "mismatch version: " << version);
}

void TrialPositionSwap::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(7934, ostr);
}

}  // namespace feasst
