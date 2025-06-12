#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "system/include/potential.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/ref_potential.h"

namespace feasst {

RefPotential::RefPotential(argtype * args) {
  class_name_ = "RefPotential";
  if (used("reference_index", *args)) {
    WARN("Deprecated RefPotential::reference_index->ref.");
  }
  reference_index_ = integer("reference_index", args, 0);
  ref_ = str("ref", args, "");
  args_ = *args;
  args->clear();
}
RefPotential::RefPotential(argtype args) : RefPotential(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(RefPotential,);

RefPotential::RefPotential(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 4017 && version <= 4018, "mismatch version: " << version);
  feasst_deserialize(&reference_index_, istr);
  if (version < 4018) {
    WARN("This restart version may not be compatible.");
    int tmp;
    feasst_deserialize(&tmp, istr);
  }
  feasst_deserialize(&args_, istr);
}

void RefPotential::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(4018, ostr);
  feasst_serialize(reference_index_, ostr);
  feasst_serialize(args_, ostr);
}

void RefPotential::run(MonteCarlo * mc) {
  mc->add_to_reference(MakePotential(args_), reference_index_, ref_);
}

}  // namespace feasst
