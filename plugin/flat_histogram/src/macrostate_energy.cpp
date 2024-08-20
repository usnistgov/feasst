#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "math/include/histogram.h"
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/criteria.h"
#include "flat_histogram/include/macrostate_energy.h"

namespace feasst {

MacrostateEnergy::MacrostateEnergy(argtype * args) :
    MacrostateEnergy(Histogram(args), args) {}
MacrostateEnergy::MacrostateEnergy(const Histogram& histogram,
    argtype * args) : Macrostate(histogram, args) {
  class_name_ = "MacrostateEnergy";
}
  //WARN("MacrostateEnergy has not been tested.");
  //num_ = ConstrainEnergy(
  //  {{"type", args_.key("particle_type").dflt("-1").str()}});
  //ASSERT(num_.type() >= -1, "particle_type: " << num_.type());
MacrostateEnergy::MacrostateEnergy(const Histogram& histogram,
    argtype args) : MacrostateEnergy(histogram, &args) {
  feasst_check_all_used(args);
}

MacrostateEnergy::MacrostateEnergy(argtype args) :
    MacrostateEnergy(Histogram(&args), &args) {
  feasst_check_all_used(args);
}

double MacrostateEnergy::value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const {
  //return num_.num_particles(system, acceptance);
  return criteria.current_energy() + acceptance.energy_new()
                                   - acceptance.energy_old();
}

class MapMacrostateEnergy {
 public:
  MapMacrostateEnergy() {
    auto hist = MakeHistogram({{"width", "1"}, {"max", "1"}});
    MacrostateEnergy(*hist).deserialize_map()["MacrostateEnergy"] =
      MakeMacrostateEnergy(*hist);
  }
};

static MapMacrostateEnergy mapper_ = MapMacrostateEnergy();

std::shared_ptr<Macrostate> MacrostateEnergy::create(std::istream& istr) const {
  return std::make_shared<MacrostateEnergy>(istr);
}

std::shared_ptr<Macrostate> MacrostateEnergy::create(argtype * args) const {
  return std::make_shared<MacrostateEnergy>(args);
}

MacrostateEnergy::MacrostateEnergy(std::istream& istr) : Macrostate(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7045, "version mismatch: " << version);
}

void MacrostateEnergy::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_macrostate_(ostr);
  feasst_serialize_version(7045, ostr);
}

}  // namespace feasst
