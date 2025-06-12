#include <cmath>
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "configuration/include/configuration.h"
#include "system/include/thermo_params.h"
#include "system/include/system.h"
#include "monte_carlo/include/rosenbluth.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/perturb.h"
#include "monte_carlo/include/trial_stage.h"

namespace feasst {

TrialStage::TrialStage(argtype * args) {
  rosenbluth_ = std::make_shared<Rosenbluth>();
  rosenbluth_->resize(integer("num_steps", args, 1));
  reference_ = integer("reference_index", args, -1);
  ref_ = str("ref", args, "");
  is_new_only_ = boolean("new_only", args, false);
}

argtype get_stage_args(argtype * args) {
  argtype tmp_args;
  for (const std::string key : {"num_steps", "reference_index", "ref", "new_only"}) {
    if (used(key, *args)) tmp_args.insert({key, str(key, args)});
  }
  return tmp_args;
}

void TrialStage::precompute(System * system) {
  select_->precompute(system);
  perturb_->precompute(select_.get(), system);
  if (!ref_.empty()) {
    const int conf = select_->configuration_index();
    reference_ = system->reference_index(conf, ref_);
  }
}

void TrialStage::before_select() {
  select_->before_select();
  perturb_->before_select();
}

bool TrialStage::select(System * system,
    Acceptance * acceptance,
    Random * random,
    TrialSelect * previous_select) {
  const int conf = select_->configuration_index();
  const bool is_selected = select_->select(acceptance->perturbed(conf),
                                           system, random, previous_select);
  DEBUG("is_selected " << is_selected);
  if (is_selected) {
    acceptance->add_to_perturbed(select_->mobile(), conf);
//    set_mobile_physical(false, system);
    DEBUG("select: " << select_->mobile().str());
    DEBUG("perturbed: " << acceptance->perturbed(conf).str());
  } else {
    acceptance->set_reject(true);
  }
  return is_selected;
}

void TrialStage::set_mobile_physical(const bool physical, System * system) {
  DEBUG("setting mobile physical " << physical);
  DEBUG("configuration_index " << select_->configuration_index());
  DEBUG("mobile " << select_->mobile().str());
  select_->get_configuration(system)->set_selection_physical(
    select_->mobile(),
    physical);
}

void TrialStage::set_rosenbluth_energy_(const int step, System * system) {
  DEBUG("select " << select_->mobile().str());
  double energy;
  if (reference_ == -1) {
    energy = system->perturbed_energy(select_->mobile(), select_->configuration_index());
  } else {
    energy = system->reference_energy(select_->mobile(), reference_, select_->configuration_index());
  }
  const double excluded = select().exclude_energy();
  DEBUG("excluded " << excluded);
  ASSERT(!std::isinf(energy), "energy: " << energy << " is inf.");
  ASSERT(!std::isnan(energy), "energy: " << energy << " is nan.");
  ASSERT(!std::isinf(excluded), "excluded: " << excluded << " is inf.");
  ASSERT(!std::isnan(excluded), "excluded: " << excluded << " is nan.");
  rosenbluth_->set_energy(step, energy, excluded);
  const int config = select_->configuration_index();
  rosenbluth_->set_energy_profile(step, system->stored_energy_profile(config));
}

void TrialStage::attempt(System * system,
    Acceptance * acceptance,
    Criteria * criteria,
    const int old,
    Random * random) {
  ASSERT(perturb_, "perturb not set");
  set_mobile_physical(true, system);
  DEBUG("setting mobile physical: " << select_->mobile().str());

  if (rosenbluth_->num() == 1) {
    select_->zero_exclude_energy();
    perturb_->perturb(system, select_.get(), random, old, acceptance);
    set_rosenbluth_energy_(0, system);
    rosenbluth_->compute(system->thermo_params().beta(), random, old);
  } else {
    for (int step = 0; step < rosenbluth_->num(); ++step) {
      // DEBUG(perturb_->class_name());
      bool is_position_held = false;
      if (step == 0 && old == 1) is_position_held = true;
      select_->zero_exclude_energy();
      perturb_->perturb(system, select_.get(), random, is_position_held, acceptance);
      DEBUG("updating state " << select_->mobile().trial_state());
      rosenbluth_->store(step, select_->mobile());
      DEBUG("ref " << reference_);
      set_rosenbluth_energy_(step, system);
      perturb_->revert(system);
    }
    rosenbluth_->compute(system->thermo_params().beta(), random, old);
    DEBUG("old " << old << " num " << rosenbluth_->num());
    if (old != 1) {
//      if (select_->is_ghost()) {
//        system->get_configuration()->revive(rosenbluth_->stored(0));
//      }
      if (rosenbluth_->chosen_step() != -1) {
        DEBUG("chosen " << rosenbluth_->chosen_step());
        DEBUG("updating positions " << rosenbluth_->chosen().str());
        DEBUG("pos0 " << rosenbluth_->chosen().site_positions()[0][0].str());
        // DEBUG("pos1 " << rosenbluth_->chosen().site_positions()[0][1].str());
        select_->get_configuration(system)->update_positions(rosenbluth_->chosen());
        // if select->is_ghost() then revive particle
      } else if (is_new_only()) {
        ASSERT(rosenbluth_->num() == 1, "assumes 1 step for mayer");
        FATAL("this should never happen");
        // system->get_configuration()->update_positions(rosenbluth_->stored(0));
      }
    }
  }
}

void TrialStage::begin_stage() {
  perturb_->begin_stage(*select_);
}

void TrialStage::mid_stage(System * system) {
  DEBUG("** mid stage **");
  select_->mid_stage();
  perturb_->mid_stage(*select_, *system);
  set_mobile_physical(false, system);
}

bool TrialStage::are_constraints_satisfied(const int old,
    const System& system) const {
  return select_->are_constraints_satisfied(old, system);
}

void TrialStage::serialize(std::ostream& ostr) const {
  feasst_serialize_version(136, ostr);
  feasst_serialize(reference_, ostr);
  feasst_serialize(ref_, ostr);
  feasst_serialize_fstdr(perturb_, ostr);
  feasst_serialize_fstdr(select_, ostr);
  feasst_serialize(rosenbluth_, ostr);
  feasst_serialize(is_new_only_, ostr);
}

TrialStage::TrialStage(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 135 && version <= 136, "version mismatch: " << version);
  feasst_deserialize(&reference_, istr);
  if (version >= 136) {
    feasst_deserialize(&ref_, istr);
  }
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
//  feasst_deserialize(rosenbluth_, istr);
// HWH for unknown reasons, this function template does not work.
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      rosenbluth_ = std::make_shared<Rosenbluth>(istr);
    }
  }
  feasst_deserialize(&is_new_only_, istr);
}

void TrialStage::set(std::shared_ptr<Perturb> perturb) { perturb_ = perturb; }

const Perturb& TrialStage::perturb() const {
  return const_cast<Perturb&>(*perturb_);
}

void TrialStage::revert(System * system) { perturb_->revert(system); }

void TrialStage::finalize(System * system) { perturb_->finalize(system); }

void TrialStage::tune(const double acceptance) { perturb_->tune(acceptance); }

std::string TrialStage::status_header() const {
  return perturb_->status_header();
}

std::string TrialStage::status() const { return perturb_->status(); }

void TrialStage::set_tunable(const double tunable) {
  perturb_->set_tunable(tunable);
}

void TrialStage::set(std::shared_ptr<TrialSelect> select) { select_ = select; }

const TrialSelect& TrialStage::select() const {
  return const_cast<TrialSelect&>(*select_);
}

TrialSelect * TrialStage::get_trial_select() { return select_.get(); }

const Rosenbluth& TrialStage::rosenbluth() const { return *rosenbluth_; }

const int TrialStage::num_steps() const { return rosenbluth_->num(); }
}  // namespace feasst
