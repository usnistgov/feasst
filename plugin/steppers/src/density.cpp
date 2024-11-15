#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/density.h"

namespace feasst {

FEASST_MAPPER(Density,);

Density::Density(argtype * args) : Analyze(args) {}
Density::Density(argtype args) : Density(&args) {
  feasst_check_all_used(args);
}

void Density::initialize(MonteCarlo * mc) {
  printer(header(*mc), output_file(mc->criteria()));
}

std::string Density::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << accumulator_->status_header() << std::endl;
  return ss.str();
}

void Density::update(const MonteCarlo& mc) {
  const System& system = mc.system();
  const double volume = configuration(system).domain().volume();
  DEBUG("volume: " << volume);
  const int num_particles = configuration(system).num_particles();
  DEBUG("num_particles: " << num_particles);
  DEBUG("state: " << state());
  accumulator_->accumulate(static_cast<double>(num_particles)/volume);
}

std::string Density::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(mc);
  }
  ss << accumulator_->status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void Density::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(8965, ostr);
}

Density::Density(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8965, "mismatch version:" << version);
}

}  // namespace feasst
