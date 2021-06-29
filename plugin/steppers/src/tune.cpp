#include "steppers/include/tune.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapTune {
 public:
  MapTune() {
    Tune().deserialize_map()["Tune"] = MakeTune();
  }
};

static MapTune mapper_ = MapTune();

Tune::Tune(argtype * args) : ModifyUpdateOnly(args) {}
Tune::Tune(argtype args) : Tune(&args) { check_all_used(args); }

void Tune::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(256, ostr);
}

Tune::Tune(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 256, "version mismatch:" << version);
}
}  // namespace feasst
