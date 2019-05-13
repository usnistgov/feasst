#include "system/include/potential.h"

namespace feasst {

void Potential::serialize(std::ostream& ostr) const {
  feasst_serialize_version(432, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize_fstdr(visit_model_, ostr);
  feasst_serialize_fstdr(model_, ostr);
  feasst_serialize(stored_energy_, ostr);
  feasst_serialize(model_params_override_, ostr);
  if (model_params_override_) {
    feasst_serialize_fstobj(model_params_, ostr);
  }
}

Potential::Potential(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(432 == version, "version mismatch: " << version);
  feasst_deserialize(&group_index_, istr);
  // feasst_deserialize_fstdr(visit_model_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      visit_model_ = visit_model_->deserialize(istr);
    }
  }
  // feasst_deserialize_fstdr(model_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      model_ = model_->deserialize(istr);
    }
  }
  feasst_deserialize(&stored_energy_, istr);
  feasst_deserialize(&model_params_override_, istr);
  if (model_params_override_) {
    feasst_deserialize_fstobj(&model_params_, istr);
  }
}

}  // namespace feasst
