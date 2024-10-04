#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "monte_carlo/include/trial_select.h"
#include "chain/include/perturb_pivot.h"

namespace feasst {

FEASST_MAPPER(PerturbPivot,);

PerturbPivot::PerturbPivot(argtype args) : PerturbPivot(&args) {
  feasst_check_all_used(args);
}
PerturbPivot::PerturbPivot(argtype * args) : PerturbRotate(args) {
  class_name_ = "PerturbPivot";
}

std::shared_ptr<Perturb> PerturbPivot::create(std::istream& istr) const {
  return std::make_shared<PerturbPivot>(istr);
}

PerturbPivot::PerturbPivot(std::istream& istr)
  : PerturbRotate(istr) {
  ASSERT(class_name_ == "PerturbPivot", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(445 == version, "mismatch version: " << version);
}

void PerturbPivot::serialize_perturb_pivot_(std::ostream& ostr) const {
  serialize_perturb_rotate_(ostr);
  feasst_serialize_version(445, ostr);
}

void PerturbPivot::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_pivot_(ostr);
}

void PerturbPivot::move(const bool is_position_held,
                        System * system,
                        TrialSelect * select,
                        Random * random,
                        Acceptance * acceptance) {
  if (is_position_held) return;
  const Position& pivot = select->anchor_position(0, 0, *system);
  DEBUG("piv " << pivot.str());
  PerturbRotate::move(system, select, random, pivot);
  DEBUG(select->mobile().site_positions()[0][0].str());
}
}  // namespace feasst
