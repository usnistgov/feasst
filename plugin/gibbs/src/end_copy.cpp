#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"
#include "gibbs/include/end_copy.h"

namespace feasst {

FEASST_MAPPER(EndCopy,);

EndCopy::EndCopy(argtype * args) {
  class_name_ = "EndCopy";
}
EndCopy::EndCopy(argtype args) : EndCopy(&args) {
  feasst_check_all_used(args);
}

EndCopy::EndCopy(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 2346 && version <= 2346, "mismatch version: " << version);
}

void EndCopy::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(2346, ostr);
}

void EndCopy::run(MonteCarlo * mc) {
  mc->set_parse_for_num_configs(1);
  mc->set_parse_replace();
}

}  // namespace feasst
