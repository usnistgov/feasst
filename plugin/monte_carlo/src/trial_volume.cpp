#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_volume.h"
#include "monte_carlo/include/trial_compute_volume.h"
#include "monte_carlo/include/trial_volume.h"

namespace feasst {

TrialVolume::TrialVolume(const argtype& args) : Trial(args) {
  add_stage(
    MakeTrialSelectParticle(args),
    MakePerturbVolume(args),
    args);
  set(MakeTrialComputeVolume());
  class_name_ = "TrialVolume";
}

class MapTrialVolume {
 public:
  MapTrialVolume() {
    auto obj = MakeTrialVolume();
    obj->deserialize_map()["TrialVolume"] = obj;
  }
};

static MapTrialVolume mapper_ = MapTrialVolume();

std::shared_ptr<Trial> TrialVolume::create(std::istream& istr) const {
  return std::make_shared<TrialVolume>(istr);
}

TrialVolume::TrialVolume(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialVolume", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(6346 == version, "mismatch version: " << version);
}

void TrialVolume::serialize_trial_volume_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(6346, ostr);
}

void TrialVolume::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_volume_(ostr);
}

}  // namespace feasst
