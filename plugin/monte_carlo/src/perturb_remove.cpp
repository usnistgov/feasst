#include "utils/include/serialize.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "monte_carlo/include/perturb_remove.h"

namespace feasst {

PerturbRemove::PerturbRemove(std::shared_ptr<Perturb> perturb) {
  class_name_ = "PerturbRemove";
  move_ = perturb;
  disable_tunable_();
}
PerturbRemove::PerturbRemove() : PerturbRemove(MakePerturbAnywhere()) {}

class MapPerturbRemove {
 public:
  MapPerturbRemove() {
    auto obj = MakePerturbRemove();
    obj->deserialize_map()["PerturbRemove"] = obj;
  }
};

static MapPerturbRemove mapper_ = MapPerturbRemove();

void PerturbRemove::perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held,
    Acceptance * acceptance) {
  DEBUG("setting finalize possible: " << select->mobile().str());
  ASSERT(select->mobile().num_particles() > 0, "error");
  set_finalize_possible(true, select);

  if (is_position_held) {
    move_->set_revert_possible(false, NULL);
  } else {
    move_->perturb(system, select, random, is_position_held, acceptance);
    move_->set_revert_possible(true, select);
  }

  // setting trial state should go last so other perturbs do not overwrite
  select->set_trial_state(2);
}

void PerturbRemove::before_select() {
  Perturb::before_select();
  move_->before_select();
}

void PerturbRemove::finalize(System * system) {
  if (finalize_possible()) {
    // HWH removals are performed by System::finalize for trial_state == 2.
    //system->finalize(finalize_select()->mobile());
  }
}

void PerturbRemove::revert(System * system) {
  DEBUG("move_->revert_possible() " << move_->revert_possible());
  move_->revert(system);
}

std::shared_ptr<Perturb> PerturbRemove::create(std::istream& istr) const {
  return std::make_shared<PerturbRemove>(istr);
}

PerturbRemove::PerturbRemove(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbRemove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(143 == version, "mismatch version: " << version);
  // HWH for unknown reasons, this function template does not work.
  // feasst_deserialize_fstdr(move_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      move_ = move_->deserialize(istr);
    }
  }
}

void PerturbRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(143, ostr);
  feasst_serialize_fstdr(move_, ostr);
}

}  // namespace feasst
