
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/analyze_factory.h"
#include "monte_carlo/include/modify_factory.h"
#include "steppers/include/seek_analyze.h"
#include "steppers/include/seek_modify.h"
#include "example/include/action_example.h"

namespace feasst {

FEASST_MAPPER(ActionExample,);

ActionExample::ActionExample(argtype * args) {
  class_name_ = "ActionExample";
  analyze_name_ = str("analyze_name", args, "");
  modify_name_ = str("modify_name", args, "");
}

ActionExample::ActionExample(argtype args) : ActionExample(&args) {
  feasst_check_all_used(args);
}

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
    ASSERT(index[0] != -1, "analyze_name:" << analyze_name_ <<
      " was not found. Is it a modify?");
    mc->get_analyze_factory()->get_analyze(index[0])->write_to_file(*mc);
  }
  if (!modify_name_.empty()) {
    const std::vector<int> index = SeekModify().index(modify_name_, *mc);
    ASSERT(index[0] != -1, "modify_name:" << modify_name_ <<
      " was not found. Is it an analyze?");
    mc->get_modify_factory()->get_modify(index[0])->write_to_file(mc);
  }
}

}  // namespace feasst
