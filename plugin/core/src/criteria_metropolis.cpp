
#include "core/include/criteria_metropolis.h"

namespace feasst {

bool CriteriaMetropolis::is_accepted(const AcceptanceCriteria accept_criteria) {
  if ( (accept_criteria.force_rejection != 1) &&
       (random_.uniform() < exp(accept_criteria.ln_metropolis_prob)) ) {
    set_current_energy(accept_criteria.energy_new);
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
