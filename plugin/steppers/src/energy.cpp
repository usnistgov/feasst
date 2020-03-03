#include "steppers/include/energy.h"

namespace feasst {

class MapEnergy {
 public:
  MapEnergy() {
    auto obj = MakeEnergy();
    obj->deserialize_map()["Energy"] = obj;
  }
};

static MapEnergy mapper_ = MapEnergy();

Energy::Energy(const argtype &args) : Analyze(args) {
  args_.init(args);
  energy_.set_block(args_.key("num_blocks").dflt(str(1e5)).integer());
}

void Energy::update(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  DEBUG("en: " << criteria->current_energy());
  DEBUG("state: " << state());
  energy_.accumulate(criteria->current_energy());
}

std::string Energy::write(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  ss << energy_.str() << " ";
  DEBUG(ss.str());
  return ss.str();
}

void Energy::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(325, ostr);
  feasst_serialize_fstobj(energy_, ostr);
}

Energy::Energy(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 325, "mismatch version:" << version);
  feasst_deserialize_fstobj(&energy_, istr);
}
}  // namespace feasst
