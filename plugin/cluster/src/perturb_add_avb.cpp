#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_select.h"
#include "cluster/include/perturb_add_avb.h"

namespace feasst {

PerturbAddAVB::PerturbAddAVB(argtype args) : PerturbAddAVB(&args) {
  feasst_check_all_used(args);
}
PerturbAddAVB::PerturbAddAVB(argtype * args) : Perturb(args) {
  class_name_ = "PerturbAddAVB";
  delay_add_ = boolean("delay_add", args, true);
  move_ = std::make_shared<PerturbMoveAVB>(args);
}

class MapPerturbAddAVB {
 public:
  MapPerturbAddAVB() {
    auto obj = MakePerturbAddAVB();
    obj->deserialize_map()["PerturbAddAVB"] = obj;
  }
};

static MapPerturbAddAVB mapper_ = MapPerturbAddAVB();

std::shared_ptr<Perturb> PerturbAddAVB::create(std::istream& istr) const {
  return std::make_shared<PerturbAddAVB>(istr);
}

void PerturbAddAVB::perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held,
    Acceptance * acceptance) {
  DEBUG("is_position_held " << is_position_held);
  //DEBUG(select->mobile().str());
  if (!delay_add_) {
    system->get_configuration()->revive(select->mobile());
  }
  move_->move(is_position_held, system, select, random, acceptance);
  move_->set_revert_possible(true, select);
  set_revert_possible(true, select);
  set_finalize_possible(true, select);
  // setting trial state should go last so other perturbs do not overwrite
  select->set_trial_state(3);
}

void PerturbAddAVB::revert(System * system) {
  DEBUG("revert_possible " << revert_possible());
  if (revert_possible()) {
//    DEBUG(revert_select()->mobile().str());
//    DEBUG("nump " << system->configuration().num_particles());
    //system->revert(revert_select()->mobile());
    move_->revert(system);
    if (!delay_add_) {
      system->get_configuration()->remove_particles(revert_select()->mobile());
    }
  }
}

void PerturbAddAVB::finalize(System * system) {
  DEBUG("finalize_possible " << finalize_possible());
  if (finalize_possible()) {
    //system->finalize(finalize_select()->mobile());
    if (delay_add_) {
      system->get_configuration()->revive(finalize_select()->mobile());
    }
  }
}

std::string PerturbAddAVB::status_header() const {
  std::stringstream ss;
  return ss.str();
}

std::string PerturbAddAVB::status() const {
  std::stringstream ss;
  return ss.str();
}

PerturbAddAVB::PerturbAddAVB(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbAddAVB", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(2908 == version, "mismatch version: " << version);
  feasst_deserialize(&delay_add_, istr);
  // HWH for unknown reasons, this function template does not work.
  // feasst_deserialize_fstdr(move_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
     move_ = std::make_shared<PerturbMoveAVB>(istr);
    }
  }
}

void PerturbAddAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(2908, ostr);
  feasst_serialize(delay_add_, ostr);
  feasst_serialize_fstdr(move_, ostr);
}

void PerturbAddAVB::precompute(TrialSelect * select, System * system) {
  select->set_ghost(true);
}

}  // namespace feasst
