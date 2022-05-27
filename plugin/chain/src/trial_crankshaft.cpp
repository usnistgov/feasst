#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "chain/include/select_segment.h"
#include "chain/include/perturb_crankshaft.h"
#include "chain/include/trial_crankshaft.h"

namespace feasst {

class MapTrialCrankshaft {
 public:
  MapTrialCrankshaft() {
    auto obj = MakeTrialCrankshaft();
    obj->deserialize_map()["TrialCrankshaft"] = obj;
  }
};

static MapTrialCrankshaft mapper_ = MapTrialCrankshaft();

TrialCrankshaft::TrialCrankshaft(argtype * args) :
  TrialMove(std::make_shared<SelectSegment>(args),
            std::make_shared<PerturbCrankshaft>(args),
            args) {
  class_name_ = "TrialCrankshaft";
  set_description("TrialCrankshaft");
}
TrialCrankshaft::TrialCrankshaft(argtype args) : TrialCrankshaft(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialCrankshaft::TrialCrankshaft(std::istream& istr) : TrialMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2948, "mismatch version: " << version);
}

void TrialCrankshaft::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(2948, ostr);
}

}  // namespace feasst
