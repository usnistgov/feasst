
#include "monte_carlo/include/criteria_metropolis.h"

namespace feasst {

bool CriteriaMetropolis::is_accepted(const Acceptance& acceptance,
    const System * system,
    const double uniform_random) {
  if ( (!acceptance.reject()) and
       (uniform_random < exp(acceptance.ln_metropolis_prob())) ) {
    DEBUG("accepted");
    set_current_energy(acceptance.energy_new());
    return true;
  }
  return false;
}

class MapCriteriaMetropolis {
 public:
  MapCriteriaMetropolis() {
    CriteriaMetropolis().deserialize_map()["CriteriaMetropolis"] = MakeCriteriaMetropolis();
  }
};

static MapCriteriaMetropolis mapper_ = MapCriteriaMetropolis();

CriteriaMetropolis::CriteriaMetropolis(std::istream& istr) : Criteria(istr) {
  feasst_deserialize_version(istr);
}

void CriteriaMetropolis::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(1, ostr);
}

}  // namespace feasst
