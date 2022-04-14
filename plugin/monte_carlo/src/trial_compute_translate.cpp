#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_compute_translate.h"

namespace feasst {

TrialComputeTranslate::TrialComputeTranslate(argtype args) : TrialComputeTranslate(&args) {
  check_all_used(args);
}
TrialComputeTranslate::TrialComputeTranslate(argtype * args) : TrialCompute(args) {
  class_name_ = "TrialComputeTranslate";
}

class MapTrialComputeTranslate {
 public:
  MapTrialComputeTranslate() {
    auto obj = MakeTrialComputeTranslate();
    obj->deserialize_map()["TrialComputeTranslate"] = obj;
  }
};

static MapTrialComputeTranslate mapper_ = MapTrialComputeTranslate();

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
    // save the chosen 'new' position
    //new_ = system->configuration().select_particle(first_sel.mobile().particle_index(0)).site(
    //  first_sel.mobile().site_index(0, 0)).position();
    if (first_stage->rosenbluth().chosen_step() != -1) {
      new_ = first_stage->rosenbluth().chosen();
      DEBUG("new " << new_.str() << " " << new_.site_positions()[0][0].str());
      // midstage will set new position as anchor
      for (TrialStage * stage : *stages) stage->mid_stage(system);
      // move selection to original position
      system->get_configuration()->update_positions(first_sel.mobile_original());
      // then, compute rosenbluth of old
      compute_rosenbluth(1, criteria, system, acceptance, stages, random);
      // finally, move back to new chosen position
      system->get_configuration()->update_positions(new_);
    }
  }
  DEBUG("New");
  DEBUG("current en: " << criteria->current_energy());
  DEBUG("old en: " << acceptance->energy_old());
  DEBUG("new en: " << acceptance->energy_new());
  DEBUG("energy change: " << acceptance->energy_new() - acceptance->energy_old());
  if ((*stages)[0]->is_new_only()) {
    acceptance->set_energy_new(acceptance->energy_new());
  } else {
    const double delta_energy = acceptance->energy_new() - acceptance->energy_old();
    acceptance->set_energy_new(criteria->current_energy() + delta_energy);
  }
}

std::shared_ptr<TrialCompute> TrialComputeTranslate::create(std::istream& istr) const {
  return std::make_shared<TrialComputeTranslate>(istr);
}

TrialComputeTranslate::TrialComputeTranslate(std::istream& istr)
  : TrialCompute(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(8958 == version, "mismatch version: " << version);
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
