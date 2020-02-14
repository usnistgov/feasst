#include "monte_carlo/include/trial_stage.h"

namespace feasst {

TrialStage::TrialStage(const argtype& args) {
  Arguments args_(args);
  args_.dont_check();
  rosenbluth_.resize(args_.key("num_steps").dflt("1").integer());
  reference_ = args_.key("reference_index").dflt("-1").integer();
  set_mayer();
}

void TrialStage::precompute(System * system) {
  select_->precompute(system);
  perturb_->precompute(select_.get(), system);
}

void TrialStage::before_select() {
  select_->before_select();
  perturb_->before_select();
}

void TrialStage::select(System * system,
    Acceptance * acceptance,
    Random * random) {
  const bool is_selected = select_->select(acceptance->perturbed(),
                                           system, random);
  if (is_selected) {
    acceptance->add_to_perturbed(select_->mobile());
    set_mobile_physical(false, system);
  } else {
    acceptance->set_reject(true);
  }
}

void TrialStage::set_mobile_physical(const bool physical, System * system) {
  system->get_configuration()->set_selection_physical(
    select_->mobile(),
    physical);
}

void TrialStage::attempt(System * system, Criteria * criteria, const int old,
    Random * random) {
  ASSERT(perturb_, "perturb not set");
  set_mobile_physical(true, system);
  for (int step = 0; step < rosenbluth_.num(); ++step) {
    // DEBUG(perturb_->class_name());
    bool is_position_held = false;
    if (step == 0 && old == 1) is_position_held = true;
    perturb_->perturb(system, select_.get(), random, is_position_held);
    rosenbluth_.store(step, select_->mobile(), system);
    if (reference_ == -1) {
      DEBUG("select " << select_->mobile().str());
      rosenbluth_.set_energy(step, system->perturbed_energy(select_->mobile()));
    } else {
      rosenbluth_.set_energy(step,
        system->reference_energy(select_->mobile(), reference_));
    }
    if (rosenbluth_.num() > 1) {
      perturb_->revert(system);
    }
  }
  rosenbluth_.compute(criteria->beta(), random);
  if (old != 1 && rosenbluth_.num() > 1) {
    if (rosenbluth_.chosen_step() != -1) {
      DEBUG("updating positions " << rosenbluth_.chosen().str());
      DEBUG("pos0 " << rosenbluth_.chosen().site_positions()[0][0].str());
      // DEBUG("pos1 " << rosenbluth_.chosen().site_positions()[0][1].str());
      system->get_configuration()->update_positions(rosenbluth_.chosen());
    } else if (is_mayer()) {
      ASSERT(rosenbluth_.num() == 1, "assumes 1 step for mayer");
      system->get_configuration()->update_positions(rosenbluth_.stored(0));
    }
  }
}

void TrialStage::mid_stage(System * system) {
  select_->mid_stage();
  set_mobile_physical(false, system);
}

bool TrialStage::are_constraints_satisfied(const System& system) const {
  return select_->are_constraints_satisfied(system);
}

void TrialStage::serialize(std::ostream& ostr) const {
  feasst_serialize_version(135, ostr);
  feasst_serialize(reference_, ostr);
  feasst_serialize_fstdr(perturb_, ostr);
  feasst_serialize_fstdr(select_, ostr);
  feasst_serialize_fstobj(rosenbluth_, ostr);
  feasst_serialize(is_mayer_, ostr);
}

TrialStage::TrialStage(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 135, "version");
  feasst_deserialize(&reference_, istr);
  // HWH for unknown reasons, this function template doesn't work
  //feasst_deserialize_fstdr(perturb_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      perturb_ = perturb_->deserialize(istr);
    }
  }
  // HWH for unknown reasons, this function template doesn't work
  //feasst_deserialize_fstdr(select_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      select_ = select_->deserialize(istr);
    }
  }
  feasst_deserialize_fstobj(&rosenbluth_, istr);
  feasst_deserialize(&is_mayer_, istr);
}

}  // namespace feasst
