#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_all.h"
#include "monte_carlo/include/perturb_volume.h"
#include "gibbs/include/compute_gibbs_volume_transfer.h"
#include "gibbs/include/trial_gibbs_volume_transfer.h"

namespace feasst {

FEASST_MAPPER(TrialGibbsVolumeTransfer,);

TrialGibbsVolumeTransfer::TrialGibbsVolumeTransfer(argtype * args) : Trial(args) {
  class_name_ = "TrialGibbsVolumeTransfer";
  set_description("TrialGibbsVolumeTransfer");
  const int to_configuration_index = integer("configuration_index0", args, 0);
  const int configuration_index = integer("configuration_index1", args, 1);
  ASSERT(configuration_index != to_configuration_index, "configuration_index0:"
    << configuration_index << " cannot equal configuration_index1:" <<
    to_configuration_index);
  const bool uniform_volume = boolean("uniform_volume", args, true);
  ASSERT(uniform_volume,
    "Gibbs volume transfers must be chosen uniformly in V.");
  args->insert({"uniform_volume", str(uniform_volume)});
  argtype args2 = *args;
  args2.insert({"configuration_index", str(to_configuration_index)});
  add_stage(
    std::make_shared<TrialSelectAll>(&args2),
    std::make_shared<PerturbVolume>(&args2),
    &args2);
  feasst_check_all_used(args2);
  args->insert({"configuration_index", str(configuration_index)});
  args->insert({"constrain_volume_change", "true"});
  add_stage(
    std::make_shared<TrialSelectAll>(args),
    std::make_shared<PerturbVolume>(args),
    args);
  set(std::make_shared<ComputeGibbsVolumeTransfer>());
}
TrialGibbsVolumeTransfer::TrialGibbsVolumeTransfer(argtype args) : TrialGibbsVolumeTransfer(&args) {
  feasst_check_all_used(args);
}

TrialGibbsVolumeTransfer::TrialGibbsVolumeTransfer(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9234, "mismatch version: " << version);
}

void TrialGibbsVolumeTransfer::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(9234, ostr);
}

}  // namespace feasst
