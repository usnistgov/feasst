#include "steppers/include/cpu_time.h"

namespace feasst {

class MapCPUTime {
 public:
  MapCPUTime() {
    CPUTime().deserialize_map()["CPUTime"] = MakeCPUTime();
  }
};

static MapCPUTime mapper_ = MapCPUTime();

CPUTime::CPUTime(const argtype& args) : AnalyzeWriteOnly(args) {}

}  // namespace feasst
