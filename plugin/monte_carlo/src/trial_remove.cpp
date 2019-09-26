#include "monte_carlo/include/trial_remove.h"

namespace feasst {

class MapTrialRemove {
 public:
  MapTrialRemove() {
    auto obj = MakeTrialRemove();
    obj->deserialize_map()["TrialRemove"] = obj;
  }
};

static MapTrialRemove mapper_ = MapTrialRemove();

std::shared_ptr<Trial> TrialRemove::create(std::istream& istr) const {
  return std::make_shared<TrialRemove>(istr);
}

TrialRemove::TrialRemove(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialRemove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(848 == version, "mismatch version: " << version);
}

void TrialRemove::serialize_trial_remove_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(848, ostr);
}

void TrialRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_remove_(ostr);
}

}  // namespace feasst
