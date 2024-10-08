#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "chain/include/trial_grow_linear.h"
#include "monte_carlo/include/perturb_distance.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

FEASST_MAPPER(TrialGrowLinear, MakeTrialComputeMove());

std::shared_ptr<Trial> TrialGrowLinear::create(std::istream& istr) const {
  return std::make_shared<TrialGrowLinear>(istr);
}

TrialGrowLinear::TrialGrowLinear(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialGrowLinear", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(469 == version, "mismatch version: " << version);
}

void TrialGrowLinear::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(469, ostr);
}

TrialGrowLinear::TrialGrowLinear(std::shared_ptr<TrialCompute> compute,
    argtype args) : Trial(&args) {
  class_name_ = "TrialGrowLinear";
  stored_args_ = args;
  set(compute);
}

void TrialGrowLinear::precompute(Criteria * criteria, System * system) {
  argtype tmp_args = stored_args_;
  const int type = integer("particle_type", &tmp_args, 0);
  tmp_args.insert({"particle_type", str(type)});
  const int num_sites = system->configuration().particle_type(type).num_sites();

  // put the first site anywhere
  argtype first_args = tmp_args;
  first_args.insert({"site", "0"});
  add_stage(
    std::make_shared<TrialSelectParticle>(&first_args),
    std::make_shared<PerturbAnywhere>(),
    &first_args);
  feasst_check_all_used(first_args);

  // for the rest, grow based on bond length only
  for (int site = 1; site < num_sites; ++site) {
    argtype args = tmp_args;
    args.insert(std::pair<std::string, std::string>("mobile_site", str(site)));
    args.insert(std::pair<std::string, std::string>("anchor_site", str(site - 1)));
    add_stage(
      std::make_shared<TrialSelectBond>(&args),
      std::make_shared<PerturbDistance>(&args),
      &args
    );
    feasst_check_all_used(args);
  }

  // precompute stages
  Trial::precompute(criteria, system);
}

}  // namespace feasst
