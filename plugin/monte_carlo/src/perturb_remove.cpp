#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_remove.h"

namespace feasst {

PerturbRemove::PerturbRemove() {
  class_name_ = "PerturbRemove";
  disable_tunable_();
}

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
    const bool is_position_held) {
  DEBUG("setting finalize possible: " << select->mobile().str());
  ASSERT(select->mobile().num_particles() > 0, "error");
  set_finalize_possible(true, select);

  if (is_position_held) {
    anywhere_.set_revert_possible(false, NULL);
  } else {
    anywhere_.perturb(system, select, random, is_position_held);
    anywhere_.set_revert_possible(true, select);
  }

  // setting trial state should go last so other perturbs do not overwrite
  select->set_trial_state(2);
}

void PerturbRemove::finalize(System * system) {
  if (finalize_possible()) {
    system->finalize(finalize_select()->mobile());
  }
}

void PerturbRemove::revert(System * system) {
  DEBUG("anywhere_.revert_possible() " << anywhere_.revert_possible());
  if (anywhere_.revert_possible()) {
    anywhere_.revert(system);
  }
}

std::shared_ptr<Perturb> PerturbRemove::create(std::istream& istr) const {
  return std::make_shared<PerturbRemove>(istr);
}

PerturbRemove::PerturbRemove(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbRemove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(143 == version, "mismatch version: " << version);
}

void PerturbRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(143, ostr);
}

}  // namespace feasst
