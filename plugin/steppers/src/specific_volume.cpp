#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/specific_volume.h"

namespace feasst {

FEASST_MAPPER(SpecificVolume,);

SpecificVolume::SpecificVolume(argtype * args) : Analyze(args) {}
SpecificVolume::SpecificVolume(argtype args) : SpecificVolume(&args) {
  feasst_check_all_used(args);
}

void SpecificVolume::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  printer(header(*mc), output_file(mc->criteria()));
}

std::string SpecificVolume::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << accumulator_->status_header() << std::endl;
  return ss.str();
}

void SpecificVolume::update(const MonteCarlo& mc) {
  const double volume = configuration(mc.system()).domain().volume();
  DEBUG("volume: " << volume);
  const int num = configuration(mc.system()).num_particles();
  DEBUG("num: " << num);
  DEBUG("state: " << state());
  if (num != 0) {
    get_accumulator()->accumulate(volume/static_cast<double>(num));
  }
}

std::string SpecificVolume::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(mc);
  }
  ss << accumulator().status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void SpecificVolume::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(8965, ostr);
}

SpecificVolume::SpecificVolume(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8965, "mismatch version:" << version);
}

}  // namespace feasst
