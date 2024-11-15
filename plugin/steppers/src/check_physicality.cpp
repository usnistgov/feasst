#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/check_physicality.h"

namespace feasst {

FEASST_MAPPER(CheckPhysicality,);

CheckPhysicality::CheckPhysicality(argtype * args) : AnalyzeUpdateOnly(args) {}
CheckPhysicality::CheckPhysicality(argtype args) : CheckPhysicality(&args) {
  feasst_check_all_used(args);
}

void CheckPhysicality::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(204, ostr);
}

CheckPhysicality::CheckPhysicality(std::istream& istr)
  : AnalyzeUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 204, "version mismatch: " << version);
}

void CheckPhysicality::update(const MonteCarlo& mc) {
  const System& system = mc.system();
  for (int ic = 0; ic < system.num_configurations(); ++ic) {
    ASSERT(system.configuration(ic).are_all_sites_physical(),
      "all sites are not physical");
  }
}
}  // namespace feasst
