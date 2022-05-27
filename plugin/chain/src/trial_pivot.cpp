#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "chain/include/select_end_segment.h"
#include "chain/include/perturb_pivot.h"
#include "chain/include/trial_pivot.h"

namespace feasst {

class MapTrialPivot {
 public:
  MapTrialPivot() {
    auto obj = MakeTrialPivot();
    obj->deserialize_map()["TrialPivot"] = obj;
  }
};

static MapTrialPivot mapper_ = MapTrialPivot();

TrialPivot::TrialPivot(argtype * args) :
  TrialMove(std::make_shared<SelectEndSegment>(args),
            std::make_shared<PerturbPivot>(args),
            args) {
  class_name_ = "TrialPivot";
  set_description("TrialPivot");
}
TrialPivot::TrialPivot(argtype args) : TrialPivot(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialPivot::TrialPivot(std::istream& istr) : TrialMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1243, "mismatch version: " << version);
}

void TrialPivot::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(1243, ostr);
}

}  // namespace feasst
