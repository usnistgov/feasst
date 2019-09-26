#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

class MapTrialComputeMove {
 public:
  MapTrialComputeMove() {
    auto obj = MakeTrialComputeMove();
    obj->deserialize_map()["TrialComputeMove"] = obj;
  }
};

static MapTrialComputeMove mapper_ = MapTrialComputeMove();

std::shared_ptr<TrialCompute> TrialComputeMove::create(std::istream& istr) const {
  return std::make_shared<TrialComputeMove>(istr);
}

TrialComputeMove::TrialComputeMove(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "TrialComputeMove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(888 == version, "mismatch version: " << version);
}

void TrialComputeMove::serialize_trial_compute_move_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(888, ostr);
}

void TrialComputeMove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_move_(ostr);
}

}  // namespace feasst
