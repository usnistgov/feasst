#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/position.h"
#include "configuration/include/select.h"
#include "configuration/include/particle.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "chain/include/end_to_end_distance.h"

namespace feasst {

FEASST_MAPPER(EndToEndDistance,);

EndToEndDistance::EndToEndDistance(argtype * args) : Analyze(args) {
  group_index_ = integer("group_index", args, 0);
}
EndToEndDistance::EndToEndDistance(argtype args) : EndToEndDistance(&args) {
  feasst_check_all_used(args);
}

void EndToEndDistance::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  printer(header(*criteria, *system, *trial_factory),
          output_file(*criteria));
}

std::string EndToEndDistance::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << accumulator().status_header() << std::endl;
  return ss.str();
}

void EndToEndDistance::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  const Select& selection = system.configuration().group_select(group_index_);
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = system.configuration().select_particle(part_index);
    const Position& site0pos = part.site(0).position();
    const Position& sitenpos = part.site(part.num_sites() - 1).position();
    const double distance = site0pos.distance(sitenpos);
    get_accumulator()->accumulate(distance);
  }
}

std::string EndToEndDistance::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(criteria, system, trial_factory);
  }
  ss << accumulator().status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void EndToEndDistance::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(7685, ostr);
  feasst_serialize(group_index_, ostr);
}

EndToEndDistance::EndToEndDistance(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7685, "mismatch version:" << version);
  feasst_deserialize(&group_index_, istr);
}

EndToEndDistance::EndToEndDistance(const Analyze& energy) {
  std::stringstream ss;
  energy.serialize(ss);
  *this = EndToEndDistance(ss);
}

}  // namespace feasst
