#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/volume.h"

namespace feasst {

FEASST_MAPPER(Volume,);

Volume::Volume(argtype * args) : Analyze(args) {}
Volume::Volume(argtype args) : Volume(&args) {
  feasst_check_all_used(args);
}

void Volume::initialize(MonteCarlo * mc) {
  printer(header(*mc), output_file(mc->criteria()));
}

std::string Volume::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << accumulator_->status_header() << std::endl;
  return ss.str();
}

void Volume::update(const MonteCarlo& mc) {
  const double volume = configuration(mc.system()).domain().volume();
  DEBUG("volume: " << volume);
  DEBUG("state: " << state());
  get_accumulator()->accumulate(volume);
}

std::string Volume::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(mc);
  }
  ss << accumulator().status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void Volume::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(8965, ostr);
}

Volume::Volume(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8965, "mismatch version:" << version);
}

}  // namespace feasst
