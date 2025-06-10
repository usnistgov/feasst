#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "system/include/potential.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/constrain_volume_by_cutoff.h"

namespace feasst {

FEASST_MAPPER(ConstrainVolumeByCutoff,);

ConstrainVolumeByCutoff::ConstrainVolumeByCutoff(argtype * args) : Constraint() {
  class_name_ = "ConstrainVolumeByCutoff";
}
ConstrainVolumeByCutoff::ConstrainVolumeByCutoff(argtype args) : ConstrainVolumeByCutoff(&args) {
  feasst_check_all_used(args);
}

bool ConstrainVolumeByCutoff::is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) {
  bool allowed = true;
  for (int iconfig = 0; iconfig < system.num_configurations(); ++iconfig) {
    const Configuration& config = system.configuration(iconfig);
    const PotentialFactory& unopt = system.unoptimized(iconfig);
    for (int ipot = 0; ipot < unopt.num(); ++ipot) {
      const Potential& pot = unopt.potential(ipot);
      if (!pot.does_cutoff_fit_domain(config)) {
        allowed = false;
      }
    }
  }
  DEBUG("allowed: " << allowed);
  return allowed;
}

ConstrainVolumeByCutoff::ConstrainVolumeByCutoff(std::istream& istr)
  : Constraint(istr) {
  // ASSERT(class_name_ == "ConstrainVolumeByCutoff", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7806 == version, "mismatch version: " << version);
}

void ConstrainVolumeByCutoff::serialize_constrain_volume_by_cutoff_(
    std::ostream& ostr) const {
  serialize_constraint_(ostr);
  feasst_serialize_version(7806, ostr);
}

void ConstrainVolumeByCutoff::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_constrain_volume_by_cutoff_(ostr);
}

}  // namespace feasst
