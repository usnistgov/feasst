#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

TrialComputeMove::TrialComputeMove(argtype args) : TrialComputeMove(&args) {
  feasst_check_all_used(args);
}
TrialComputeMove::TrialComputeMove(argtype * args) : TrialCompute(args) {
  class_name_ = "TrialComputeMove";
}

FEASST_MAPPER(TrialComputeMove,);

void TrialComputeMove::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeMove");
  for (TrialStage * stage : *stages) stage->begin_stage();
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  for (TrialStage * stage : *stages) stage->mid_stage(system);
  DEBUG("New");
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  const int config = stages->front()->select().configuration_index();
  DEBUG("config " << config);
  DEBUG("current en: " << criteria->current_energy(config));
  DEBUG("old en: " << acceptance->energy_old(config));
  DEBUG("new en: " << acceptance->energy_new(config));
  DEBUG("energy change: " << acceptance->energy_new(config) - acceptance->energy_old(config));
  if ((*stages)[0]->is_new_only()) {
    //acceptance->set_energy_new(acceptance->energy_new());
  } else {
    const double delta_energy = acceptance->energy_new(config) - acceptance->energy_old(config);
    acceptance->set_energy_new(criteria->current_energy(config) + delta_energy, config);
    acceptance->add_to_energy_profile_new(criteria->current_energy_profile(config), config);
    acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old(config), config);
  }
}

std::shared_ptr<TrialCompute> TrialComputeMove::create(std::istream& istr) const {
  return std::make_shared<TrialComputeMove>(istr);
}

TrialComputeMove::TrialComputeMove(std::istream& istr)
  : TrialCompute(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(888 == version, "mismatch version: " << version);
}

void TrialComputeMove::serialize_trial_compute_move_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(888, ostr);
}

void TrialComputeMove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_move_(ostr);
}

}  // namespace feasst
