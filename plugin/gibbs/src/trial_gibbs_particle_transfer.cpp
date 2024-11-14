#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "gibbs/include/compute_gibbs_particle_transfer.h"
#include "gibbs/include/trial_gibbs_particle_transfer.h"

namespace feasst {

FEASST_MAPPER(TrialGibbsParticleTransferOneWay,
  argtype({{"to_configuration_index", "1"}}));

TrialGibbsParticleTransferOneWay::TrialGibbsParticleTransferOneWay(argtype * args) : Trial(args) {
  class_name_ = "TrialGibbsParticleTransferOneWay";
  set_description("TrialGibbsParticleTransferOneWay");
  const int to_configuration_index = integer("to_configuration_index", args);
  const int configuration_index = integer("configuration_index", args, 0);
  argtype add_args = *args;
  add_args.insert({"configuration_index", str(to_configuration_index)});
  add_stage(
    std::make_shared<TrialSelectParticle>(&add_args),
    std::make_shared<PerturbAdd>(),
    &add_args);
  feasst_check_all_used(add_args);
  args->insert({"configuration_index", str(configuration_index)});
  add_stage(
    std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbRemove>(),
    args);
  set(std::make_shared<ComputeGibbsParticleTransfer>());
}
TrialGibbsParticleTransferOneWay::TrialGibbsParticleTransferOneWay(argtype args) : TrialGibbsParticleTransferOneWay(&args) {
  feasst_check_all_used(args);
}

TrialGibbsParticleTransferOneWay::TrialGibbsParticleTransferOneWay(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2409, "mismatch version: " << version);
}

void TrialGibbsParticleTransferOneWay::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(2409, ostr);
}

FEASST_MAPPER(TrialGibbsParticleTransfer,);

TrialGibbsParticleTransfer::TrialGibbsParticleTransfer(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialGibbsParticleTransfer";
  const int config0 = integer("configuration_index0", args, 0);
  const int config1 = integer("configuration_index1", args, 1);
  //INFO("args " << str(*args));
  //ASSERT(!used("configuration_index", *args),
  //  "Do not use argument:configuration_index. Use configuration_index0 or 1.");
  argtype orig_args = *args;
  orig_args.insert({"configuration_index", str(config0)});
  orig_args.insert({"to_configuration_index", str(config1)});
  args->insert({"configuration_index", str(config1)});
  args->insert({"to_configuration_index", str(config0)});
  auto trial_add = MakeTrialGibbsParticleTransferOneWay(orig_args);
  trial_add->set_weight(trial_add->weight()/2.);
  add(trial_add);
  auto trial_remove = std::make_shared<TrialGibbsParticleTransferOneWay>(args);
  trial_remove->set_weight(trial_remove->weight()/2.);
  add(trial_remove);
}
TrialGibbsParticleTransfer::TrialGibbsParticleTransfer(argtype args) : TrialGibbsParticleTransfer(&args) {
  feasst_check_all_used(args);
}

}  // namespace feasst
