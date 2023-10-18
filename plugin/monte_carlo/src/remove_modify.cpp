
#include <fstream>
#include <sstream>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "system/include/visit_model_cell.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/remove_modify.h"

namespace feasst {

RemoveModify::RemoveModify(argtype * args) {
  class_name_ = "RemoveModify";
  index_ = integer("index", args, -1);
  name_ = str("name", args, "");
  all_ = boolean("all", args, false);
}
RemoveModify::RemoveModify(argtype args) : RemoveModify(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapRemoveModify {
 public:
  MapRemoveModify() {
    auto obj = MakeRemoveModify();
    obj->deserialize_map()["RemoveModify"] = obj;
  }
};

static MapRemoveModify mapper_RemoveModify = MapRemoveModify();

RemoveModify::RemoveModify(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2045, "mismatch version: " << version);
  feasst_deserialize(&index_, istr);
  feasst_deserialize(&name_, istr);
  feasst_deserialize(&all_, istr);
}

void RemoveModify::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(2045, ostr);
  feasst_serialize(index_, ostr);
  feasst_serialize(name_, ostr);
  feasst_serialize(all_, ostr);
}

void RemoveModify::run(MonteCarlo * mc) {
  DEBUG("name " << name_);
  if (!name_.empty()) {
    for (int modify = 0; modify < mc->num_modifiers(); ++modify) {
      DEBUG("mod " << mc->modify(modify).class_name());
      if (mc->modify(modify).class_name() == name_) {
        ASSERT(index_ < 0 || modify == index_,
          "RemoveModify cannot specify both index and name");
        index_ = modify;
        DEBUG("removing " << modify);
        break;
      }
    }
    ASSERT(index_ != -1, "No Modify of name: " << name_);
  }
  DEBUG("index " << index_);
  if (index_ >= 0) {
    DEBUG("removing " << index_);
    mc->remove_modify(index_);
  }
  if (all_) {
    for (int i = 0; mc->num_modifiers() > 0; ++i) {
      mc->remove_modify(0);
      ASSERT(i < 1e5, "Infinite loop?");
    }
  }
}

}  // namespace feasst
