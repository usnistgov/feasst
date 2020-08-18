#include <cmath>
#include "utils/include/serialize.h"
#include "confinement/include/always_accept.h"

namespace feasst {

AlwaysAccept::AlwaysAccept(const argtype &args) : Criteria(args) {
  class_name_ = "AlwaysAccept";
}

AlwaysAccept::AlwaysAccept(std::shared_ptr<Constraint> constraint,
    const argtype& args) : AlwaysAccept(args) {
  add(constraint);
}

class MapAlwaysAccept {
 public:
  MapAlwaysAccept() {
    AlwaysAccept().deserialize_map()["AlwaysAccept"] = MakeAlwaysAccept();
  }
};

static MapAlwaysAccept mapper_ = MapAlwaysAccept();

AlwaysAccept::AlwaysAccept(std::istream& istr) : Criteria(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4204, "version mismatch: " << version);
}

void AlwaysAccept::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(4204, ostr);
}

bool AlwaysAccept::is_accepted(const Acceptance& acceptance,
    const System& system,
    Random * random) {
  set_current_energy(acceptance.energy_new());
  return true;
}

}  // namespace feasst
