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

WallClockLimit::WallClockLimit(argtype * args) : AnalyzeUpdateOnly(args) {
  max_hours_ = dble("max_hours", args);
  set_steps_per(integer("steps_per", args, 1));
}
WallClockLimit::WallClockLimit(argtype args) : WallClockLimit(&args) {
  check_all_used(args);
}

void WallClockLimit::update(const Criteria& criteria,
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
