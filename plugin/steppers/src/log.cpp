#include "steppers/include/log.h"
#include "utils/include/serialize.h"

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

void Log::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  std::stringstream ss;
  ss << system->status_header()
     << criteria->status_header()
     << trial_factory->status_header()
     << std::endl;
  printer(ss.str());
}

std::string Log::write(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  // ensure the following order matches the header from initialization.
  std::stringstream ss;
  ss << system.status()
     << criteria->status()
     << trial_factory.status()
     << std::endl;
  return ss.str();
}

void Log::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(668, ostr);
}

Log::Log(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 668, "version mismatch:" << version);
}

}  // namespace feasst
