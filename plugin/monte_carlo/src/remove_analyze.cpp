#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "system/include/visit_model_cell.h"
#include "monte_carlo/include/analyze_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/remove_analyze.h"

namespace feasst {

RemoveAnalyze::RemoveAnalyze(argtype * args) {
  class_name_ = "RemoveAnalyze";
  index_ = integer("index", args, -1);
  name_ = str("name", args, "");
  all_ = boolean("all", args, false);
}
RemoveAnalyze::RemoveAnalyze(argtype args) : RemoveAnalyze(&args) {
  feasst_check_all_used(args);
}

class MapRemoveAnalyze {
 public:
  MapRemoveAnalyze() {
    auto obj = MakeRemoveAnalyze();
    obj->deserialize_map()["RemoveAnalyze"] = obj;
  }
};

static MapRemoveAnalyze mapper_RemoveAnalyze = MapRemoveAnalyze();

RemoveAnalyze::RemoveAnalyze(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7985, "mismatch version: " << version);
  feasst_deserialize(&index_, istr);
  feasst_deserialize(&name_, istr);
  feasst_deserialize(&all_, istr);
}

void RemoveAnalyze::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(7985, ostr);
  feasst_serialize(index_, ostr);
  feasst_serialize(name_, ostr);
  feasst_serialize(all_, ostr);
}

void RemoveAnalyze::run(MonteCarlo * mc) {
  DEBUG("name " << name_);
  if (!name_.empty()) {
    for (int analyze = 0; analyze < mc->num_analyzers(); ++analyze) {
      DEBUG("an " << mc->analyze(analyze).class_name());
      if (mc->analyze(analyze).class_name() == name_) {
        ASSERT(index_ < 0 || analyze == index_,
          "RemoveAnalyze cannot specify both index and name");
        index_ = analyze;
        DEBUG("removing " << analyze);
        break;
      }
    }
    ASSERT(index_ != -1, "No Analyze of name: " << name_);
  }
  DEBUG("index " << index_);
  if (index_ >= 0) {
    DEBUG("removing " << index_);
    mc->remove_analyze(index_);
  }
  if (all_) {
    for (int i = 0; mc->num_analyzers() > 0; ++i) {
      mc->remove_analyze(0);
      ASSERT(i < 1e5, "Infinite loop?");
    }
  }
}

}  // namespace feasst
