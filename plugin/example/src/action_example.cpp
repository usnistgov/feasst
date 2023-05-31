
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/seek_analyze.h"
#include "example/include/action_example.h"

namespace feasst {

ActionExample::ActionExample(argtype * args) {
  class_name_ = "ActionExample";
  analyze_name_ = str("analyze_name", args, "");
  modify_name_ = str("modify_name", args, "");
}
ActionExample::ActionExample(argtype args) : ActionExample(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapActionExample {
 public:
  MapActionExample() {
    auto obj = MakeActionExample();
    obj->deserialize_map()["ActionExample"] = obj;
  }
};

static MapActionExample mapper_ActionExample = MapActionExample();

ActionExample::ActionExample(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4932, "mismatch version: " << version);
  feasst_deserialize(&analyze_name_, istr);
  feasst_deserialize(&modify_name_, istr);
}

void ActionExample::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(4932, ostr);
  feasst_serialize(analyze_name_, ostr);
  feasst_serialize(modify_name_, ostr);
}

void ActionExample::run(MonteCarlo * mc) {
  if (!analyze_name_.empty()) {
    const std::vector<int> index = SeekAnalyze().index(analyze_name_, *mc);
    ASSERT(index[1] == -1, "ActionExample not implemented for multistate");
    INFO(feasst_str(index));
  }
}

}  // namespace feasst
