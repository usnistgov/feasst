#include "monte_carlo/include/trial_rotate.h"

namespace feasst {

TrialRotate::TrialRotate(const argtype& args)
  : TrialMove(
    std::make_shared<TrialSelectParticle>(),
    std::make_shared<PerturbRotate>(args),
    args
  ) {
  class_name_ = "TrialRotate";
}

class MapTrialRotate {
 public:
  MapTrialRotate() {
    auto obj = MakeTrialRotate();
    obj->deserialize_map()["TrialRotate"] = obj;
  }
};

static MapTrialRotate mapper_ = MapTrialRotate();

std::shared_ptr<Trial> TrialRotate::create(std::istream& istr) const {
  return std::make_shared<TrialRotate>(istr);
}

TrialRotate::TrialRotate(std::istream& istr)
  : TrialMove(istr) {
  // ASSERT(class_name_ == "TrialRotate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(169 == version, "mismatch version: " << version);
}

void TrialRotate::serialize_trial_rotate_(std::ostream& ostr) const {
  serialize_trial_move_(ostr);
  feasst_serialize_version(169, ostr);
}

void TrialRotate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_rotate_(ostr);
}

}  // namespace feasst
