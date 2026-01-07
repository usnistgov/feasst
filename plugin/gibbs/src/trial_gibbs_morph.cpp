#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "gibbs/include/perturb_particle_type.h"
#include "gibbs/include/compute_gibbs_morph.h"
#include "gibbs/include/trial_gibbs_morph.h"

namespace feasst {

FEASST_MAPPER(TrialGibbsMorphOneWay, argtype({{"to_config", "1"},
  {"particle_type", "0"}, {"particle_type_morph", "1"}}));

TrialGibbsMorphOneWay::TrialGibbsMorphOneWay(argtype * args) : Trial(args) {
  class_name_ = "TrialGibbsMorphOneWay";
  set_description("TrialGibbsMorphOneWay");
  const std::string to_config = str("to_config", args);
  const std::string config = str("config", args, "0");
  DEBUG("to_config: " << to_config);
  DEBUG("config: " << config);
  DEBUG("args: " << str(*args));
  const std::string t1 = str("particle_type", args);
  const std::string t2 = str("particle_type_morph", args);
  DEBUG("t1: " << t1 << " t2: " << t2);
  argtype first_args = *args;
  first_args.insert({"config", to_config});
  first_args.insert({"particle_type", t1});
  add_stage(
    std::make_shared<TrialSelectParticle>(&first_args),
    std::make_shared<PerturbParticleType>(argtype({{"type", t2}})),
    &first_args);
  feasst_check_all_used(first_args);
  args->insert({"config", config});
  args->insert({"particle_type", t2});
  add_stage(
    std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbParticleType>(argtype({{"type", t1}})),
    args);
  set(std::make_shared<ComputeGibbsMorph>());
}
TrialGibbsMorphOneWay::TrialGibbsMorphOneWay(argtype args) : TrialGibbsMorphOneWay(&args) {
  feasst_check_all_used(args);
}

TrialGibbsMorphOneWay::TrialGibbsMorphOneWay(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4864, "mismatch version: " << version);
}

void TrialGibbsMorphOneWay::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(4864, ostr);
}

FEASST_MAPPER(TrialGibbsMorph, argtype({{"particle_type", "0"},
                                        {"particle_type_morph", "1"}}));

TrialGibbsMorph::TrialGibbsMorph(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialGibbsMorph";
  std::string cargs = str("configs", args, "0,1");
  std::vector<std::string> configs = split(cargs, ',');
  ASSERT(static_cast<int>(configs.size()) == 2, "Requires two configuration " <<
    "names separated by a comma (e.g., \"vapor,liquid\") but was given:" <<
    cargs);
  argtype orig_args = *args;
  orig_args.insert({"config", configs[0]});
  orig_args.insert({"to_config", configs[1]});
  args->insert({"config", configs[1]});
  args->insert({"to_config", configs[0]});
  auto trial_first = std::make_shared<TrialGibbsMorphOneWay>(orig_args);
  trial_first->set_weight(trial_first->weight()/2.);
  add(trial_first);
  auto trial_second = std::make_shared<TrialGibbsMorphOneWay>(args);
  trial_second->set_weight(trial_second->weight()/2.);
  add(trial_second);
}
TrialGibbsMorph::TrialGibbsMorph(argtype args) : TrialGibbsMorph(&args) {
  feasst_check_all_used(args);
}

}  // namespace feasst
