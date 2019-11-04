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
  if (args_.key("append").dflt("true").boolean()) {
    set_append();
  } else {
    ERROR("append is required");
  }
}

}  // namespace feasst
