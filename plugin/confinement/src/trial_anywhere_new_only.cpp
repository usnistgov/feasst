#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "mayer/include/trial.h"
#include "confinement/include/trial_anywhere_new_only.h"

namespace feasst {

TrialAnywhereNewOnly::TrialAnywhereNewOnly(const argtype& args)
  : TrialMove(
    std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbAnywhere>(),
    args) {
  class_name_ = "TrialAnywhereNewOnly";
  set_new_only(true);
  set(std::make_shared<TrialComputeMoveNewOnly>());
}

class MapTrialAnywhereNewOnly {
 public:
  MapTrialAnywhereNewOnly() {
    auto obj = MakeTrialAnywhereNewOnly();
    obj->deserialize_map()["TrialAnywhereNewOnly"] = obj;
  }
};

static MapTrialAnywhereNewOnly mapper_ = MapTrialAnywhereNewOnly();

std::shared_ptr<Trial> TrialAnywhereNewOnly::create(std::istream& istr) const {
  return std::make_shared<TrialAnywhereNewOnly>(istr);
}

TrialAnywhereNewOnly::TrialAnywhereNewOnly(std::istream& istr) : TrialMove(istr) {
  // ASSERT(class_name_ == "TrialAnywhereNewOnly", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(1411 == version, "mismatch version: " << version);
}

void TrialAnywhereNewOnly::serialize_trial_anywhere_new_only_(
    std::ostream& ostr) const {
  serialize_trial_move_(ostr);
  feasst_serialize_version(1411, ostr);
}

void TrialAnywhereNewOnly::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_anywhere_new_only_(ostr);
}

}  // namespace feasst
