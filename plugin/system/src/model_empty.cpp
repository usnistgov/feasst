#include "utils/include/serialize.h"
#include "system/include/model_empty.h"

namespace feasst {

FEASST_MAPPER(ModelEmpty,);

void ModelEmpty::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(189, ostr);
}

ModelEmpty::ModelEmpty(std::istream& istr) : ModelOneBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(189 == version, version);
}

double ModelEmpty::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) {
  FATAL("Empty model should not be called");
}

}  // namespace feasst
