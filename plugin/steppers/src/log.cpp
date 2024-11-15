#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/log.h"

namespace feasst {

FEASST_MAPPER(Log,);

Log::Log(argtype args) : Log(&args) { feasst_check_all_used(args); }
Log::Log(argtype * args) : AnalyzeWriteOnly(args) {
  if (boolean("append", args, true)) {
    set_append();
  } else {
    ERROR("append is required");
  }
  max_precision_ = boolean("max_precision", args, false);
  include_bonds_ = boolean("include_bonds", args, true);
}

void Log::initialize(MonteCarlo * mc) {
  printer(header(*mc), output_file(mc->criteria()));
}

std::string Log::header(const MonteCarlo& mc) const {
  const System& system = mc.system();
  std::stringstream ss;
  ss << system.status_header();
  ss << mc.criteria().status_header(system);
  if (include_bonds_) {
    std::string append = "";
    for (int iconf = 0; iconf < system.num_configurations(); ++iconf) {
      if (system.num_configurations() > 1) {
        append = "_config" + str(iconf);
      }
      ss << ",BondTwoBody" << append <<
            ",BondThreeBody" << append <<
            ",BondFourBody" << append;
    }
  }
  // print number of trials here instead of TrialFactory header because
  // multiple factories makes it redundant.
  ss << ",trial"
     << mc.trial_factory().status_header()
     << std::endl;
  return ss.str();
}

std::string Log::write(const MonteCarlo& mc) {
  const System& system = mc.system();
  const TrialFactory& trial_factory = mc.trial_factory();
  // ensure the following order matches the header from initialization.
  std::stringstream ss;
  ss << system.status();
  ss << mc.criteria().status(max_precision_);
  if (include_bonds_) {
    for (int iconf = 0; iconf < system.num_configurations(); ++iconf) {
      bond_visitor_.compute_all(system.configuration());
      if (max_precision_) {
        ss << "," << MAX_PRECISION << bond_visitor_.energy_two_body()
           << "," << MAX_PRECISION << bond_visitor_.energy_three_body()
           << "," << MAX_PRECISION << bond_visitor_.energy_four_body();
      } else {
        ss << "," << bond_visitor_.energy_two_body()
           << "," << bond_visitor_.energy_three_body()
           << "," << bond_visitor_.energy_four_body();
      }
    }
  }
  ss << "," << trial_factory.num_attempts()
     << trial_factory.status()
     << std::endl;
  return ss.str();
}

void Log::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(668, ostr);
  feasst_serialize(max_precision_, ostr);
  feasst_serialize(include_bonds_, ostr);
  feasst_serialize_fstobj(bond_visitor_, ostr);
}

Log::Log(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 668, "version mismatch:" << version);
  feasst_deserialize(&max_precision_, istr);
  feasst_deserialize(&include_bonds_, istr);
  feasst_deserialize_fstobj(&bond_visitor_, istr);
}

}  // namespace feasst
