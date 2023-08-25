#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "monte_carlo/include/metropolis.h"

namespace feasst {

Metropolis::Metropolis(argtype * args) : Criteria(args) {
  class_name_ = "Metropolis";
  // HWH depreciate
  // Support depreciation warning for old argument name
  if (used("num_attempts_per_iteration", *args)) {
    WARN("Metropolis argument num_attempts_per_iteration is depreciated. " <<
         "Use num_trials_per_iteration instead.");
    ASSERT(!used("num_trials_per_iteration", *args),
      "Both num_trials_per_iteration and num_attempts_per_iteration");
    num_trials_per_iteration_ =
      integer("num_attempts_per_iteration", args);
  }
  num_trials_per_iteration_ =
    integer("num_trials_per_iteration", args, 1e9);
}
Metropolis::Metropolis(argtype args) : Metropolis(&args) {
  FEASST_CHECK_ALL_USED(args);
}

Metropolis::Metropolis(std::shared_ptr<Constraint> constraint) : Metropolis() {
  add(constraint);
}

bool Metropolis::is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) {
  check_num_iterations_(num_trials_per_iteration_);
  DEBUG("ln_prob " << acceptance->ln_metropolis_prob());
  DEBUG("is_allowed " << is_allowed(system, *acceptance));
  DEBUG("reject " << acceptance->reject());
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
  feasst_deserialize(&num_trials_per_iteration_, istr);
}

void Metropolis::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(278, ostr);
  feasst_serialize(num_trials_per_iteration_, ostr);
}

}  // namespace feasst
