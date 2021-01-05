#include <cmath>
#include "utils/include/serialize.h"
#include "confinement/include/henry_coefficient.h"

namespace feasst {

class MapHenryCoefficient {
 public:
  MapHenryCoefficient() {
    auto obj = MakeHenryCoefficient();
    obj->deserialize_map()["HenryCoefficient"] = obj;
  }
};

static MapHenryCoefficient mapper_ = MapHenryCoefficient();

HenryCoefficient::HenryCoefficient(const argtype &args) : Analyze(args) {
  args_.init(args);
  ASSERT(steps_per_update() == 1, "should update every step");
}

void HenryCoefficient::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  ASSERT(criteria->class_name() == "AlwaysReject",
    "HenryCoefficient requires AlwaysReject");
  ASSERT(trial_factory->num() == 1,
    "HenryCoefficient requires only one Trial");
  ASSERT(trial_factory->trial(0).description() == "TrialAdd",
    "HenryCoefficient requires TrialAdd. " <<
    "Found: " << trial_factory->trial(0).description());
  ASSERT(trial_factory->trial(0).compute().is_new_only(),
    "HenryCoefficient requires new_only.");
  printer(header(*criteria, *system, *trial_factory),
          file_name(*criteria));
}

std::string HenryCoefficient::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << accumulator_.status_header() << std::endl;
  return ss.str();
}

void HenryCoefficient::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  const double en = trial_factory.trial(0).accept().energy_new();
  DEBUG("en: " << en);
  accumulator_.accumulate(std::exp(-system.thermo_params().beta()*en));
}

std::string HenryCoefficient::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  ss << accumulator_.status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void HenryCoefficient::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(9492, ostr);
}

HenryCoefficient::HenryCoefficient(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9492, "mismatch version:" << version);
}

}  // namespace feasst
