#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/histogram.h"
#include "configuration/include/configuration.h"
#include "configuration/include/particle_factory.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "flat_histogram/include/macrostate_position.h"

namespace feasst {

MacrostatePosition::MacrostatePosition(argtype * args) :
    MacrostatePosition(Histogram(args), args) {}
MacrostatePosition::MacrostatePosition(const Histogram& histogram,
    argtype * args) : Macrostate(histogram, args) {
  class_name_ = "MacrostatePosition";
  particle_index_ = integer("particle_index", args, 0);
  site_index_ = integer("site_index", args, 0);
  dimension_ = integer("dimension", args, 0);
}
MacrostatePosition::MacrostatePosition(const Histogram& histogram,
    argtype args) : MacrostatePosition(histogram, &args) {
  feasst_check_all_used(args);
}
MacrostatePosition::MacrostatePosition(argtype args) :
    MacrostatePosition(Histogram(&args), &args) {
  feasst_check_all_used(args);
}
std::shared_ptr<Macrostate> MacrostatePosition::create(argtype * args) const {
  return std::make_shared<MacrostatePosition>(args);
}

double MacrostatePosition::value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const {
  return system.configuration().particle(particle_index_).site(site_index_).position().coord(dimension_);
}

FEASST_MAPPER(MacrostatePosition, argtype({{"width", "1"}, {"max", "1"}}));

std::shared_ptr<Macrostate> MacrostatePosition::create(std::istream& istr) const {
  return std::make_shared<MacrostatePosition>(istr);
}

MacrostatePosition::MacrostatePosition(std::istream& istr)
  : Macrostate(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8274, "version mismatch: " << version);
  feasst_deserialize(&particle_index_, istr);
  feasst_deserialize(&site_index_, istr);
  feasst_deserialize(&dimension_, istr);
}

void MacrostatePosition::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_macrostate_(ostr);
  feasst_serialize_version(8274, ostr);
  feasst_serialize(particle_index_, ostr);
  feasst_serialize(site_index_, ostr);
  feasst_serialize(dimension_, ostr);
}

}  // namespace feasst
