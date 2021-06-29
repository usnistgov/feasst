#include <iostream>  // std::scientific
#include "utils/include/serialize.h"
#include "steppers/include/cpu_time.h"

namespace feasst {

class MapCPUTime {
 public:
  MapCPUTime() {
    CPUTime().deserialize_map()["CPUTime"] = MakeCPUTime();
  }
};

static MapCPUTime mapper_ = MapCPUTime();

CPUTime::CPUTime(argtype * args) : AnalyzeWriteOnly(args) {}
CPUTime::CPUTime(argtype args) : CPUTime(&args) {
  check_all_used(args);
}

void CPUTime::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  num_writes_ = 0;
  initialize_time_ = cpu_hours();
}

std::string CPUTime::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  ++num_writes_;
  const double elapsed_hours = cpu_hours() - initialize_time_;
  const double steps_per_second = num_writes_*steps_per_write()
                                  /(elapsed_hours*60*60);
  accumulator_.accumulate(steps_per_second);
  ss << std::scientific
     << "elapsed_hours: " << elapsed_hours << " "
     << "steps per second: " << steps_per_second << " "
     << std::fixed
     << std::endl;
  return ss.str();
}

void CPUTime::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(235, ostr);
  feasst_serialize(num_writes_, ostr);
  feasst_serialize(initialize_time_, ostr);
}

CPUTime::CPUTime(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(235 == version, "version mismatch:" << version);
  feasst_deserialize(&num_writes_, istr);
  feasst_deserialize(&initialize_time_, istr);
}

}  // namespace feasst
