#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "system/include/potential.h"
#include "system/include/model.h"
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"
#include "model_expanded/include/constrain_model_index.h"

namespace feasst {

ConstrainModelIndex::ConstrainModelIndex(argtype * args) : Constraint() {
  class_name_ = "ConstrainModelIndex";
  maximum_ = integer("maximum", args, -1);
  minimum_ = integer("minimum", args, 0);
  potential_index_ = integer("potential_index", args, 0);
}
ConstrainModelIndex::ConstrainModelIndex(argtype args) : ConstrainModelIndex(&args) {
  feasst_check_all_used(args);
}

int ConstrainModelIndex::model_index(const System& system,
    const Acceptance& acceptance) const {
  const int shift = acceptance.macrostate_shift();
  DEBUG("shift " << shift);
  DEBUG("potential_index_ " << potential_index_);
  return system.potential(potential_index_).model().model_index() + shift;
}

bool ConstrainModelIndex::is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const {
  const int index = model_index(system, acceptance);
  bool allowed;
  if (index >= minimum_ && ( (maximum_ == -1) || (index <= maximum_) ) ) {
    allowed = true;
  } else {
    allowed = false;
  }
  DEBUG("index: " << index << " allowed: " << allowed << " max " << maximum_ <<
    " min " << minimum_);
  return allowed;
}

FEASST_MAPPER(ConstrainModelIndex,);

ConstrainModelIndex::ConstrainModelIndex(std::istream& istr)
  : Constraint(istr) {
  // ASSERT(class_name_ == "ConstrainModelIndex", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(4503 == version, "mismatch version: " << version);
  feasst_deserialize(&maximum_, istr);
  feasst_deserialize(&minimum_, istr);
  feasst_deserialize(&potential_index_, istr);
}

void ConstrainModelIndex::serialize_constrain_model_index_(
    std::ostream& ostr) const {
  serialize_constraint_(ostr);
  feasst_serialize_version(4503, ostr);
  feasst_serialize(maximum_, ostr);
  feasst_serialize(minimum_, ostr);
  feasst_serialize(potential_index_, ostr);
}

void ConstrainModelIndex::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_constrain_model_index_(ostr);
}

}  // namespace feasst
