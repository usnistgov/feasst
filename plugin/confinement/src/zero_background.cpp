#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "system/include/potential.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "confinement/include/background.h"
#include "confinement/include/zero_background.h"

namespace feasst {

ZeroBackground::ZeroBackground(argtype * args) {
  class_name_ = "ZeroBackground";
  configuration_index_ = integer("configuration_index", args, 0);
}
ZeroBackground::ZeroBackground(argtype args) : ZeroBackground(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(ZeroBackground,);

ZeroBackground::ZeroBackground(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1497, "mismatch version: " << version);
  feasst_deserialize(&configuration_index_, istr);
}

void ZeroBackground::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(1497, ostr);
  feasst_serialize(configuration_index_, ostr);
}

void ZeroBackground::run(MonteCarlo * mc) {
  const double current_energy = mc->get_system()->energy();
  auto bg = MakeBackground({{"constant", str(-current_energy)}});
  mc->add(std::make_shared<Potential>(bg), configuration_index_);
}

}  // namespace feasst
