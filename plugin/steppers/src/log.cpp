#include "steppers/include/log.h"

namespace feasst {

class MapLog {
 public:
  MapLog() {
    Log().deserialize_map()["Log"] = MakeLog();
  }
};

static MapLog mapper_ = MapLog();

Log::Log(const argtype& args) : AnalyzeWriteOnly(args) {
  set_append();
}

}  // namespace feasst
