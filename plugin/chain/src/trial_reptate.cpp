#include "chain/include/trial_reptate.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapTrialReptate {
 public:
  MapTrialReptate() {
    auto obj = MakeTrialReptate({{"max_length", "1"}});
    obj->deserialize_map()["TrialReptate"] = obj;
  }
};

static MapTrialReptate mapper_ = MapTrialReptate();

std::shared_ptr<Trial> TrialReptate::create(std::istream& istr) const {
  return std::make_shared<TrialReptate>(istr);
}

TrialReptate::TrialReptate(std::istream& istr) : TrialMove(istr) {
  // ASSERT(class_name_ == "TrialReptate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(159 == version, "mismatch version: " << version);
}

void TrialReptate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(159, ostr);
}

}  // namespace feasst
