#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "configuration/include/neighbor_criteria.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_select.h"
#include "cluster/include/perturb_move_avb.h"

namespace feasst {

PerturbMoveAVB::PerturbMoveAVB(argtype args) : PerturbMoveAVB(&args) {
  feasst_check_all_used(args);
}
PerturbMoveAVB::PerturbMoveAVB(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbMoveAVB";
  neighbor_ = integer("neighbor_index", args, 0);
  rotate_.set_tunable(180.);
  rotate_.disable_tunable_();
  disable_tunable_();
  inside_ = boolean("inside", args, true);
}

FEASST_MAPPER(PerturbMoveAVB,);

std::shared_ptr<Perturb> PerturbMoveAVB::create(std::istream& istr) const {
  return std::make_shared<PerturbMoveAVB>(istr);
}

void PerturbMoveAVB::move(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    Acceptance * acceptance) {
  if (is_position_held) return;
  DEBUG("mobile " << select->mobile().str());
  DEBUG("anchor " << select->anchor().str());
  const Configuration& config = system->configuration();

  // give particle random orientation
  if (select->mobile().num_sites() > 1) {
    rotate_.move(is_position_held, system, select, random, acceptance);
  }

  if (displace_.dimension() == 0) {
    displace_.set_to_origin(config.dimension());
  }
  const Position& anchor_pos = config.select_particle(
    select->anchor().particle_index(0)).site(
    select->anchor().site_index(0, 0)).position();
  DEBUG("anchor_pos " << anchor_pos.str());

  DEBUG("inside " << inside_);
  NeighborCriteria * neighbor = system->get_neighbor_criteria(neighbor_, select->configuration_index());
  if (inside_) {
    // translate mobile particle to AV of anchor
    random->position_in_spherical_shell(
      neighbor->minimum_distance(),
      neighbor->maximum_distance(),
      &displace_);
  } else {
    // pick a random position in the domain but outside the AV, considering PBC
    const double volume_av = neighbor->volume(config.dimension());
    ASSERT(config.domain().volume() > volume_av, "AV: " << volume_av
      << " too large for domain");
    int attempts = 0;
    bool found = false;
    while (!found) {
      DEBUG("attempts: " << attempts);
      config.domain().random_position(&displace_, random);
      displace_.subtract(anchor_pos);
      found = !neighbor->is_position_accepted(displace_, config.domain());
      ++attempts;
      ASSERT(attempts < 1e5, "max attempts");
    }
  }
  DEBUG("displace: " << displace_.str());
  // move from anchor_pos reference to mobile reference
  displace_.add(anchor_pos);
  DEBUG("displace + anchor: " << displace_.str());
  DEBUG("mobile: " << select->mobile().str());
  const Position& mobile_pos = select->mobile().site_positions()[0][0];
  DEBUG("mobile_pos " << mobile_pos.str());
  displace_.subtract(mobile_pos);
  DEBUG("displace + anchor - mobile" << displace_.str());
  DEBUG("old pos " << config.select_particle(select->mobile().particle_index(0)).site(0).position().str());
  translate_.move(displace_, system, select);
  DEBUG("new pos " << config.select_particle(select->mobile().particle_index(0)).site(0).position().str());
}

PerturbMoveAVB::PerturbMoveAVB(std::istream& istr)
  : PerturbMove(istr) {
  ASSERT(class_name_ == "PerturbMoveAVB", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(5937 == version, "mismatch version: " << version);
  feasst_deserialize(&neighbor_, istr);
  feasst_deserialize_fstobj(&rotate_, istr);
  feasst_deserialize(&inside_, istr);
}

void PerturbMoveAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(5937, ostr);
  feasst_serialize(neighbor_, ostr);
  feasst_serialize_fstobj(rotate_, ostr);
  feasst_serialize(inside_, ostr);
}

}  // namespace feasst
