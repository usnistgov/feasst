#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "monte_carlo/include/metropolis.h"

namespace feasst {

Metropolis::Metropolis(argtype args) : Criteria(&args) {
  class_name_ = "Metropolis";
  num_attempts_per_iteration_ =
    integer("num_attempts_per_iteration", &args, 1e9);
  check_all_used(args);
}

Metropolis::Metropolis(std::shared_ptr<Constraint> constraint) : Metropolis() {
  add(constraint);
}

bool Metropolis::is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) {
  check_num_iterations_(num_attempts_per_iteration_);
  DEBUG("ln_prob " << acceptance->ln_metropolis_prob());
  if ( (!acceptance->reject()) &&
       is_allowed(system, *acceptance) &&
       (random->uniform() < std::exp(acceptance->ln_metropolis_prob())) ) {
    DEBUG("accepted");
    set_current_energy(acceptance->energy_new());
    set_current_energy_profile(acceptance->energy_profile_new());
    was_accepted_ = true;
  } else {
    DEBUG("rejected");
    was_accepted_ = false;
  }
  acceptance->set_endpoint(false);
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
  feasst_deserialize(&num_attempts_per_iteration_, istr);
}

void Metropolis::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(278, ostr);
  feasst_serialize(num_attempts_per_iteration_, ostr);
}

}  // namespace feasst
