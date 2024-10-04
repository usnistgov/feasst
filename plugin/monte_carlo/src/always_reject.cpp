#include "utils/include/serialize.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/always_reject.h"

namespace feasst {

FEASST_MAPPER(AlwaysReject,);

AlwaysReject::AlwaysReject() : Criteria() {
  class_name_ = "AlwaysReject";
}

AlwaysReject::AlwaysReject(std::shared_ptr<Constraint> constraint)
  : AlwaysReject() {
  add(constraint);
}

AlwaysReject::AlwaysReject(std::istream& istr) : Criteria(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4204, "version mismatch: " << version);
}

void AlwaysReject::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(4204, ostr);
}

bool AlwaysReject::is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) {
  DEBUG("en new " << acceptance->energy_new());
  DEBUG("en old " << acceptance->energy_old());
  DEBUG("en new prof " << feasst_str(acceptance->energy_profile_new()));
  return false;
}

}  // namespace feasst
