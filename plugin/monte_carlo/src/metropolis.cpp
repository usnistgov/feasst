
#include "monte_carlo/include/metropolis.h"

namespace feasst {

bool Metropolis::is_accepted(const Acceptance& acceptance,
    const System * system,
    const double uniform_random) {
  if ( (!acceptance.reject()) and
       (uniform_random < exp(acceptance.ln_metropolis_prob())) ) {
    DEBUG("accepted");
    set_current_energy(acceptance.energy_new());
    was_accepted_ = true;
  } else {
    DEBUG("rejected");
    was_accepted_ = false;
  }
  return was_accepted_;
}

class MapMetropolis {
 public:
  MapMetropolis() {
    Metropolis().deserialize_map()["Metropolis"] = MakeMetropolis();
  }
};

static MapMetropolis mapper_ = MapMetropolis();

Metropolis::Metropolis(std::istream& istr) : Criteria(istr) {
  feasst_deserialize_version(istr);
}

void Metropolis::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(1, ostr);
}

}  // namespace feasst
