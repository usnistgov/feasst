#include "utils/include/serialize.h"
#include "steppers/include/check.h"

namespace feasst {

class MapCheck {
 public:
  MapCheck() {
    Check().deserialize_map()["Check"] = MakeCheck();
  }
};

static MapCheck mapper_ = MapCheck();

void Check::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(890, ostr);
}

Check::Check(std::istream& istr) : AnalyzeUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 890, "unrecognized version: " << version);
}

void Check::update(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  system.configuration().check();
  system.check();
}

}  // namespace feasst
