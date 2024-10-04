#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "monte_carlo/include/trial_select.h"
#include "chain/include/perturb_particle_pivot.h"

namespace feasst {

PerturbParticlePivot::PerturbParticlePivot(argtype args) : PerturbParticlePivot(&args) {
  feasst_check_all_used(args);
}
PerturbParticlePivot::PerturbParticlePivot(argtype * args) : PerturbRotate(args) {
  class_name_ = "PerturbParticlePivot";
}

FEASST_MAPPER(PerturbParticlePivot,);

std::shared_ptr<Perturb> PerturbParticlePivot::create(std::istream& istr) const {
  return std::make_shared<PerturbParticlePivot>(istr);
}

void PerturbParticlePivot::move(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    Acceptance * acceptance) {
  if (is_position_held) return;
  ASSERT(select->mobile().num_sites() > 0, "selection error");
  ASSERT(select->mobile().site_positions().size() > 0, "requires coordinates");
  DEBUG("num particle positions " << select->mobile().site_positions().size());
  DEBUG("num site positions " << select->mobile().site_positions()[0].size());
  ASSERT(select->anchor().num_sites() == 1, "selection error");
  const Position& pivot = select->anchor_position(0, 0, *system);
  DEBUG("pivot " << pivot.str());
  PerturbRotate::move(system, select, random, pivot);
}

PerturbParticlePivot::PerturbParticlePivot(std::istream& istr)
  : PerturbRotate(istr) {
  // ASSERT(class_name_ == "PerturbParticlePivot", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(4398 == version, "mismatch version: " << version);
}

void PerturbParticlePivot::serialize_perturb_particle_pivot_(std::ostream& ostr) const {
  serialize_perturb_rotate_(ostr);
  feasst_serialize_version(4398, ostr);
}

void PerturbParticlePivot::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_particle_pivot_(ostr);
}

}  // namespace feasst
