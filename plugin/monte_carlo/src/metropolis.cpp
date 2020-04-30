#include <cmath>
#include "monte_carlo/include/metropolis.h"
#include "utils/include/serialize.h"

namespace feasst {

Metropolis::Metropolis(const argtype &args) : Criteria(args) {
  class_name_ = "Metropolis";
}

Metropolis::Metropolis(std::shared_ptr<Constraint> constraint,
    const argtype& args) : Metropolis(args) {
  add(constraint);
}

bool Metropolis::is_accepted(const Acceptance& acceptance,
    const System * system,
    const double uniform_random) {
  if ( (!acceptance.reject()) &&
       (is_allowed(system, acceptance)) &&
       (uniform_random < std::exp(acceptance.ln_metropolis_prob())) ) {
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
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 278, "version mismatch: " << version);
}

void Metropolis::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(278, ostr);
}

}  // namespace feasst
