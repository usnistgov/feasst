#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "system/include/potential.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/opt_potential.h"

namespace feasst {

OptPotential::OptPotential(argtype * args) {
  class_name_ = "OptPotential";
  configuration_index_ = integer("configuration_index", args, 0);
  args_ = *args;
  args->clear();
}
OptPotential::OptPotential(argtype args) : OptPotential(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(OptPotential,);

OptPotential::OptPotential(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6070, "mismatch version: " << version);
  feasst_deserialize(&configuration_index_, istr);
  feasst_deserialize(&args_, istr);
}

void OptPotential::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(6070, ostr);
  feasst_serialize(configuration_index_, ostr);
  feasst_serialize(args_, ostr);
}

void OptPotential::run(MonteCarlo * mc) {
  mc->add_to_optimized(MakePotential(args_),
                       configuration_index_);
}

}  // namespace feasst
