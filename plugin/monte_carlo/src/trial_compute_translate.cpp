#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/rosenbluth.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_compute_translate.h"

namespace feasst {

TrialComputeTranslate::TrialComputeTranslate(argtype args) : TrialComputeTranslate(&args) {
  feasst_check_all_used(args);
}
TrialComputeTranslate::TrialComputeTranslate(argtype * args) : TrialCompute(args) {
  class_name_ = "TrialComputeTranslate";
  new_ = std::make_shared<Select>();
}

FEASST_MAPPER(TrialComputeTranslate,);

void TrialComputeTranslate::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeTranslate");
  TrialStage * first_stage = (*stages)[0];
  //const Perturb& first_pert = first_stage->perturb();
  const TrialSelect& first_sel = first_stage->select();
  //ASSERT(first_pert.class_name() == "PerturbTranslate", "error");

  for (TrialStage * stage : *stages) stage->begin_stage();
  if (first_stage->num_steps() == 1) {
    compute_rosenbluth(1, criteria, system, acceptance, stages, random);
    for (TrialStage * stage : *stages) stage->mid_stage(system);
    compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  } else {
    //ASSERT(first_sel.mobile().num_sites() == 1, "multi stage translate only " <<
    //  "implemented for 1 site");
    DEBUG("original pos " << first_sel.mobile_original().site_positions()[0][0].str());
    // compute rosenbluth of new first
    compute_rosenbluth(0, criteria, system, acceptance, stages, random);
    if (first_stage->rosenbluth().chosen_step() != -1) {
      new_ = std::make_shared<Select>(first_stage->rosenbluth().chosen());
      DEBUG("new " << new_->str() << " " << new_->site_positions()[0][0].str());
      // midstage will set new position as anchor
      for (TrialStage * stage : *stages) stage->mid_stage(system);
      // move selection to original position
      system->get_configuration()->update_positions(first_sel.mobile_original());
      // then, compute rosenbluth of old
      compute_rosenbluth(1, criteria, system, acceptance, stages, random);
      // finally, move back to new chosen position
      system->get_configuration()->update_positions(*new_);
    }
  }
  if (!acceptance->reject()) {
    ASSERT(system->num_configurations() == 1, "not implemented for multiple configs");
    DEBUG("New");
    DEBUG("current en: " << criteria->current_energy());
    DEBUG("old en: " << acceptance->energy_old());
    DEBUG("new en: " << acceptance->energy_new());
    DEBUG("energy change: " << acceptance->energy_new() - acceptance->energy_old());
    if ((*stages)[0]->is_new_only()) {
      //acceptance->set_energy_new(acceptance->energy_new());
    } else {
      const double delta_energy = acceptance->energy_new() - acceptance->energy_old();
      acceptance->set_energy_new(criteria->current_energy() + delta_energy);
      acceptance->add_to_energy_profile_new(criteria->current_energy_profile());
      acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old());
    }
  }
}

std::shared_ptr<TrialCompute> TrialComputeTranslate::create(std::istream& istr) const {
  return std::make_shared<TrialComputeTranslate>(istr);
}

TrialComputeTranslate::TrialComputeTranslate(std::istream& istr)
  : TrialCompute(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(8958 == version, "mismatch version: " << version);
  new_ = std::shared_ptr<Select>();
}

void TrialComputeTranslate::serialize_trial_compute_translate_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(8958, ostr);
}

void TrialComputeTranslate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_translate_(ostr);
}

}  // namespace feasst
