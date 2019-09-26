#include "monte_carlo/include/trial_compute_add.h"

namespace feasst {

class MapTrialComputeAdd {
 public:
  MapTrialComputeAdd() {
    auto obj = MakeTrialComputeAdd();
    obj->deserialize_map()["TrialComputeAdd"] = obj;
  }
};

static MapTrialComputeAdd mapper_ = MapTrialComputeAdd();

std::shared_ptr<TrialCompute> TrialComputeAdd::create(std::istream& istr) const {
  return std::make_shared<TrialComputeAdd>(istr);
}

TrialComputeAdd::TrialComputeAdd(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "TrialComputeAdd", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(834 == version, "mismatch version: " << version);
}

void TrialComputeAdd::serialize_trial_compute_add_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(834, ostr);
}

void TrialComputeAdd::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_add_(ostr);
}

}  // namespace feasst
