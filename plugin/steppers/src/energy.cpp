#include "steppers/include/energy.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapEnergy {
 public:
  MapEnergy() {
    auto obj = MakeEnergy();
    obj->deserialize_map()["Energy"] = obj;
  }
};

static MapEnergy mapper_ = MapEnergy();

Energy::Energy(argtype args) : Analyze(&args) {
  check_all_used(args);
}

void Energy::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  printer(header(*criteria, *system, *trial_factory),
          file_name(*criteria));
}

std::string Energy::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << accumulator_.status_header() << std::endl;
  return ss.str();
}

void Energy::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  DEBUG("en: " << criteria.current_energy());
  DEBUG("state: " << state());
  accumulator_.accumulate(criteria.current_energy());
}

std::string Energy::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  ss << accumulator_.status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void Energy::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(325, ostr);
}

Energy::Energy(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 325, "mismatch version:" << version);
}

Energy::Energy(const Analyze& energy) {
  std::stringstream ss;
  energy.serialize(ss);
  *this = Energy(ss);
}

}  // namespace feasst
