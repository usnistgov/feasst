#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/metropolis.h"

namespace feasst {

FEASST_MAPPER(Metropolis,);

Metropolis::Metropolis(argtype * args) : Criteria(args) {
  class_name_ = "Metropolis";
  // HWH deprecate
  // Support deprecation warning for old argument name
  if (used("num_trials_per_iteration", *args)) {
    WARN("Metropolis argument num_trials_per_iteration is deprecated. " <<
         "Use trials_per_cycle instead.");
    trials_per_cycle_ = integer("num_trials_per_iteration", args);
  } else {
    trials_per_cycle_ = integer("trials_per_cycle", args, 1e9);
  }
  DEBUG("trials_per_cycle " << trials_per_cycle_);
}
Metropolis::Metropolis(argtype args) : Metropolis(&args) {
  feasst_check_all_used(args);
}

Metropolis::Metropolis(std::shared_ptr<Constraint> constraint) : Metropolis() {
  add(constraint);
}

bool Metropolis::is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) {
  check_num_cycles_(trials_per_cycle_);
  DEBUG("is_allowed " << is_allowed(system, *acceptance));
  DEBUG("reject " << acceptance->reject());
  DEBUG("ln_prob " << acceptance->ln_metropolis_prob());
  if ( (!acceptance->reject()) &&
       is_allowed(system, *acceptance) &&
       (random->uniform() < std::exp(acceptance->ln_metropolis_prob())) ) {
    DEBUG("accepted");
    update_current_energy(*acceptance);
    was_accepted_ = true;
  } else {
    DEBUG("rejected");
    was_accepted_ = false;
  }
  acceptance->set_endpoint(false);
  return was_accepted_;
}

Metropolis::Metropolis(std::istream& istr) : Criteria(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 278, "version mismatch: " << version);
  feasst_deserialize(&trials_per_cycle_, istr);
}

void Metropolis::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(278, ostr);
  feasst_serialize(trials_per_cycle_, ostr);
}

}  // namespace feasst
