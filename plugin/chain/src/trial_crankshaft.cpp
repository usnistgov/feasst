#include "chain/include/trial_crankshaft.h"
#include "utils/include/serialize.h"
#include "chain/include/select_segment.h"
#include "chain/include/perturb_crankshaft.h"

namespace feasst {

class MapTrialCrankshaft {
 public:
  MapTrialCrankshaft() {
    auto obj = MakeTrialCrankshaft();
    obj->deserialize_map()["TrialCrankshaft"] = obj;
  }
};

static MapTrialCrankshaft mapper_ = MapTrialCrankshaft();

TrialCrankshaft::TrialCrankshaft(const argtype& args) : TrialMove(
    std::make_shared<SelectSegment>(args),
    std::make_shared<PerturbCrankshaft>(args),
    args
  ) {
  class_name_ = "TrialCrankshaft";
}

std::shared_ptr<Trial> TrialCrankshaft::create(std::istream& istr) const {
  return std::make_shared<TrialCrankshaft>(istr);
}

TrialCrankshaft::TrialCrankshaft(std::istream& istr) : TrialMove(istr) {
  // ASSERT(class_name_ == "TrialCrankshaft", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(121 == version, "mismatch version: " << version);
}

void TrialCrankshaft::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(121, ostr);
}

}  // namespace feasst
