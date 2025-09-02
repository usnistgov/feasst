#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/specific_energy.h"

namespace feasst {

FEASST_MAPPER(SpecificEnergy,);

SpecificEnergy::SpecificEnergy(argtype * args) : Analyze(args) {}
SpecificEnergy::SpecificEnergy(argtype args) : SpecificEnergy(&args) {
  feasst_check_all_used(args);
}

void SpecificEnergy::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  printer(header(*mc), output_file(mc->criteria()));
}

std::string SpecificEnergy::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << accumulator_->status_header() << std::endl;
  return ss.str();
}

void SpecificEnergy::update(const MonteCarlo& mc) {
  const double en = mc.criteria().current_energy(configuration_index());
  DEBUG("en: " << en);
  const int num = configuration(mc.system()).num_particles();
  DEBUG("num: " << num);
  DEBUG("state: " << state());
  if (num != 0) {
    accumulator_->accumulate(en/static_cast<double>(num));
  }
}

std::string SpecificEnergy::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(mc);
  }
  ss << accumulator_->status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void SpecificEnergy::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(325, ostr);
}

SpecificEnergy::SpecificEnergy(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 325, "mismatch version:" << version);
}

}  // namespace feasst
