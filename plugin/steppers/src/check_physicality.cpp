#include "steppers/include/check_physicality.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapCheckPhysicality {
 public:
  MapCheckPhysicality() {
    CheckPhysicality().deserialize_map()["CheckPhysicality"] =
      MakeCheckPhysicality();
  }
};

static MapCheckPhysicality mapper_ = MapCheckPhysicality();

CheckPhysicality::CheckPhysicality(argtype * args) : AnalyzeUpdateOnly(args) {}
CheckPhysicality::CheckPhysicality(argtype args) : CheckPhysicality(&args) {
  FEASST_CHECK_ALL_USED(args);
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

void CheckPhysicality::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  for (int ic = 0; ic < system.num_configurations(); ++ic) {
    ASSERT(system.configuration(ic).are_all_sites_physical(),
      "all sites are not physical");
  }
}
}  // namespace feasst
