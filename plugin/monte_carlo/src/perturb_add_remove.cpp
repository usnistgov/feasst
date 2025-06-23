#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/perturb_add_remove.h"

namespace feasst {

FEASST_MAPPER(PerturbAddRemove,);

PerturbAddRemove::PerturbAddRemove(argtype * args) : Perturb(args) {
  class_name_ = "PerturbAddRemove";
  add_ = std::make_unique<PerturbAdd>(args);
  rm_ = std::make_unique<PerturbRemove>();
}
PerturbAddRemove::PerturbAddRemove(argtype args) : PerturbAddRemove(&args) {
  feasst_check_all_used(args);
}
PerturbAddRemove::~PerturbAddRemove() {}

std::shared_ptr<Perturb> PerturbAddRemove::create(std::istream& istr) const {
  return std::make_shared<PerturbAddRemove>(istr);
}

void PerturbAddRemove::precompute(TrialSelect * select, System * system) {
  add_->precompute(select, system);
  rm_->precompute(select, system);
}

void PerturbAddRemove::before_select() {
  Perturb::before_select();
  add_->before_select();
  rm_->before_select();
}

void PerturbAddRemove::perturb(System * system, TrialSelect * select, Random * random,
    const bool is_position_held, Acceptance * acceptance) {
  set_finalize_possible(true, select);
  Perturb * perturb;
  if (select->is_ghost()) {
    adding_ = true;
    perturb = add_.get();
  } else {
    adding_ = false;
    perturb = rm_.get();
  }
  DEBUG("adding:" << adding_);
  perturb->perturb(system, select, random, is_position_held, acceptance);
  if (is_position_held) {
    set_revert_possible(false, NULL);
  } else {
    set_revert_possible(true, select);
  }
}

void PerturbAddRemove::revert(System * system) {
  if (revert_possible()) {
    if (adding_) {
      add_->revert(system);
    } else {
      rm_->revert(system);
    }
  }
}

void PerturbAddRemove::finalize(System * system) {
  if (finalize_possible()) {
    if (adding_) {
      add_->finalize(system);
    } else {
      rm_->finalize(system);
    }
  }
}

std::string PerturbAddRemove::status_header() const {
  std::stringstream ss;
  return ss.str();
}

std::string PerturbAddRemove::status() const {
  std::stringstream ss;
  return ss.str();
}

PerturbAddRemove::PerturbAddRemove(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbAddRemove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7439 == version, "mismatch version: " << version);
  feasst_deserialize(add_, istr);
  feasst_deserialize(rm_, istr);
}

void PerturbAddRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(7439, ostr);
  feasst_serialize(add_, ostr);
  feasst_serialize(rm_, ostr);
}

}  // namespace feasst
