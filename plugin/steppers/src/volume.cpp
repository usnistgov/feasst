#include "configuration/include/domain.h"
#include "steppers/include/volume.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapVolume {
 public:
  MapVolume() {
    auto obj = MakeVolume();
    obj->deserialize_map()["Volume"] = obj;
  }
};

static MapVolume mapper_ = MapVolume();

Volume::Volume(argtype * args) : Analyze(args) {}
Volume::Volume(argtype args) : Volume(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void Volume::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  printer(header(*criteria, *system, *trial_factory),
          file_name(*criteria));
}

std::string Volume::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << accumulator_.status_header() << std::endl;
  return ss.str();
}

void Volume::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  const double volume = system.configuration(configuration()).domain().volume();
  DEBUG("volume: " << volume);
  DEBUG("state: " << state());
  accumulator_.accumulate(volume);
}

std::string Volume::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  ss << accumulator_.status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void Volume::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(8965, ostr);
}

Volume::Volume(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8965, "mismatch version:" << version);
}

}  // namespace feasst
