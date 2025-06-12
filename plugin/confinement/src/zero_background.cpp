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
  if (used("configuration_index", *args)) {
    WARN("Deprecated ZeroBackground::configuration_index->config (see Configuration::name)");
  }
  configuration_index_ = integer("configuration_index", args, 0);
  config_ = str("config", args, "");
}
ZeroBackground::ZeroBackground(argtype args) : ZeroBackground(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(ZeroBackground,);

ZeroBackground::ZeroBackground(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 1497 && version <= 1498, "mismatch version: " << version);
  feasst_deserialize(&configuration_index_, istr);
  if (version >= 1498) {
    feasst_deserialize(&config_, istr);
  }
}

void ZeroBackground::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(1498, ostr);
  feasst_serialize(configuration_index_, ostr);
  feasst_serialize(config_, ostr);
}

void ZeroBackground::run(MonteCarlo * mc) {
  if (!config_.empty()) {
    configuration_index_ = mc->system().configuration_index(config_);
  }
  const double current_energy = mc->get_system()->energy();
  auto bg = MakeBackground({{"constant", str(-current_energy)}});
  mc->add(std::make_shared<Potential>(bg,
    argtype({{"configuration_index", str(configuration_index_)},
             {"config", config_}})));
}

}  // namespace feasst
