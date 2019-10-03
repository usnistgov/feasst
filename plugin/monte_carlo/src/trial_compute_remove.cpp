#include "monte_carlo/include/trial_compute_remove.h"

namespace feasst {

TrialComputeRemove::TrialComputeRemove() {
  class_name_ = "TrialComputeRemove";
}

class MapTrialComputeRemove {
 public:
  MapTrialComputeRemove() {
    auto obj = MakeTrialComputeRemove();
    obj->deserialize_map()["TrialComputeRemove"] = obj;
  }
};

static MapTrialComputeRemove mapper_ = MapTrialComputeRemove();

std::shared_ptr<TrialCompute> TrialComputeRemove::create(std::istream& istr) const {
  return std::make_shared<TrialComputeRemove>(istr);
}

TrialComputeRemove::TrialComputeRemove(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "TrialComputeRemove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(332 == version, "mismatch version: " << version);
}

void TrialComputeRemove::serialize_trial_compute_remove_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(332, ostr);
}

void TrialComputeRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_remove_(ostr);
}

}  // namespace feasst
