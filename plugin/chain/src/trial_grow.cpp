#include "utils/include/serialize.h"
#include "chain/include/trial_grow.h"
#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

class MapTrialGrowLinear {
 public:
  MapTrialGrowLinear() {
    auto obj = MakeTrialGrowLinear(MakeTrialComputeMove());
    obj->deserialize_map()["TrialGrowLinear"] = obj;
  }
};

static MapTrialGrowLinear mapper_ = MapTrialGrowLinear();

std::shared_ptr<Trial> TrialGrowLinear::create(std::istream& istr) const {
  return std::make_shared<TrialGrowLinear>(istr);
}

TrialGrowLinear::TrialGrowLinear(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialGrowLinear", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(469 == version, "mismatch version: " << version);
}

void TrialGrowLinear::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(469, ostr);
}

}  // namespace feasst
