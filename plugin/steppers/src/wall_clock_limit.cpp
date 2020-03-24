#include "steppers/include/wall_clock_limit.h"
#include "utils/include/serialize.h"

namespace feasst {

// this example shows how to handle a required argument.
class MapWallClockLimit {
 public:
  MapWallClockLimit() {
    auto obj = MakeWallClockLimit({{"max_hours", "0"}});
    obj->deserialize_map()["WallClockLimit"] = obj;
  }
};

static MapWallClockLimit mapper_ = MapWallClockLimit();

WallClockLimit::WallClockLimit(const argtype &args) : AnalyzeUpdateOnly(args) {
  args_.init(args);
  max_hours_ = args_.key("max_hours").dble();
  set_steps_per(args_.key("steps_per").dflt("1").integer());
}

void WallClockLimit::update(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  const double hours = double(clock())/double(CLOCKS_PER_SEC)/60./60.;
  ASSERT(hours < max_hours_, "wall clock hours(" << hours << ") exceed " <<
    "the maximum(" << max_hours_ << ")");
}

void WallClockLimit::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(895, ostr);
  feasst_serialize(max_hours_, ostr);
}

WallClockLimit::WallClockLimit(std::istream& istr) : AnalyzeUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 895, "version mismatch:" << version);
  feasst_deserialize(&max_hours_, istr);
}

}  // namespace feasst
