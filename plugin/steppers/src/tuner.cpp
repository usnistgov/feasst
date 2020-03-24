#include "steppers/include/tuner.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapTuner {
 public:
  MapTuner() {
    Tuner().deserialize_map()["Tuner"] = MakeTuner();
  }
};

static MapTuner mapper_ = MapTuner();

void Tuner::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(256, ostr);
}

Tuner::Tuner(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 256, "version mismatch:" << version);
}
}  // namespace feasst
