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
  argtype({{"to_config", "1"}}));

TrialGibbsParticleTransferOneWay::TrialGibbsParticleTransferOneWay(argtype * args) : Trial(args) {
  class_name_ = "TrialGibbsParticleTransferOneWay";
  set_description("TrialGibbsParticleTransferOneWay");
  const std::string to_config = str("to_config", args);
  const std::string config = str("config", args, "0");
  argtype add_args = *args;
  add_args.insert({"config", to_config});
  add_stage(
    std::make_shared<TrialSelectParticle>(&add_args),
    std::make_shared<PerturbAdd>(),
    &add_args);
  feasst_check_all_used(add_args);
  args->insert({"config", str(config)});
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
  std::string cargs = str("configs", args, "0,1");
  std::vector<std::string> configs = split(cargs, ',');
  ASSERT(static_cast<int>(configs.size()) == 2, "Requires two configuration " <<
    "names separated by a comma (e.g., \"vapor,liquid\") but was given:" <<
    cargs);
  //INFO("args " << str(*args));
  //ASSERT(!used("configuration_index", *args),
  //  "Do not use argument:configuration_index. Use configuration_index0 or 1.");
  argtype orig_args = *args;
  orig_args.insert({"config", configs[0]});
  orig_args.insert({"to_config", configs[1]});
  args->insert({"config", configs[1]});
  args->insert({"to_config", configs[0]});
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
