#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/energy.h"

namespace feasst {

FEASST_MAPPER(Energy,);

Energy::Energy(argtype * args) : Analyze(args) {}
Energy::Energy(argtype args) : Energy(&args) {
  feasst_check_all_used(args);
}

void Energy::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  printer(header(*mc), output_file(mc->criteria()));
}

std::string Energy::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << accumulator_->status_header() << std::endl;
  return ss.str();
}

void Energy::update(const MonteCarlo& mc) {
  const double en = mc.criteria().current_energy(configuration_index());
  DEBUG("en: " << en);
  DEBUG("state: " << state());
  accumulator_->accumulate(en);
}

std::string Energy::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(mc);
  }
  ss << accumulator_->status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void Energy::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(325, ostr);
}

Energy::Energy(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 325, "mismatch version:" << version);
}

Energy::Energy(const Analyze& energy) {
  std::stringstream ss;
  energy.serialize(ss);
  *this = Energy(ss);
}

}  // namespace feasst
