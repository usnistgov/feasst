
#include <fstream>
#include <sstream>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "system/include/visit_model_cell.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/remove_trial.h"

namespace feasst {

RemoveTrial::RemoveTrial(argtype * args) {
  class_name_ = "RemoveTrial";
  index_ = integer("index", args, -1);
  name_ = str("name", args, "");
  all_ = boolean("all", args, false);
  name_contains_ = str("name_contains", args, "");
}
RemoveTrial::RemoveTrial(argtype args) : RemoveTrial(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapRemoveTrial {
 public:
  MapRemoveTrial() {
    auto obj = MakeRemoveTrial();
    obj->deserialize_map()["RemoveTrial"] = obj;
  }
};

static MapRemoveTrial mapper_RemoveTrial = MapRemoveTrial();

RemoveTrial::RemoveTrial(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 3854 && version <= 3855, "mismatch version: " << version);
  feasst_deserialize(&index_, istr);
  feasst_deserialize(&name_, istr);
  if (version >= 3855) {
    feasst_deserialize(&name_contains_, istr);
  }
  feasst_deserialize(&all_, istr);
}

void RemoveTrial::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3855, ostr);
  feasst_serialize(index_, ostr);
  feasst_serialize(name_, ostr);
  feasst_serialize(name_contains_, ostr);
  feasst_serialize(all_, ostr);
}

void RemoveTrial::run(MonteCarlo * mc) {
  if (!name_.empty()) {
    for (int trial = 0; trial < mc->trials().num(); ++trial) {
      if (mc->trial(trial).class_name() == name_) {
        ASSERT(index_ < 0 || trial == index_,
          "RemoveTrial cannot specify both index and name");
        index_ = trial;
        break;
      }
    }
    ASSERT(index_ != -1, "No Trial of name: " << name_);
  }
  if (index_ >= 0) {
    mc->remove_trial(index_);
  }
  if (!name_contains_.empty()) {
    for (int trial = mc->trials().num() - 1; trial >= 0; --trial) {
      std::string name = mc->trial(trial).class_name();
      if (name == "Trial") {
        name = mc->trial(trial).description();
      }
      if (name.find(name_contains_) != std::string::npos) {
        mc->remove_trial(trial);
      }
    }
  }
  if (all_) {
    for (int i = 0; mc->trials().num() > 0; ++i) {
      mc->remove_trial(0);
      ASSERT(i < 1e5, "too many trials. Infinite loop?");
    }
  }
}

}  // namespace feasst
