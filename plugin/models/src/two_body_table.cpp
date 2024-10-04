#include "utils/include/serialize.h"
#include "models/include/two_body_table.h"

namespace feasst {

FEASST_MAPPER(TwoBodyTable,);

void TwoBodyTable::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_ideal_gas_(ostr);
}

double TwoBodyTable::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  FATAL("Should be used in conjunction with TablePotential " <<
    "(e.g., Potential Model TwoBodyTable VisitModelInner TablePotential");
}

void TwoBodyTable::serialize_ideal_gas_(std::ostream& ostr) const {
  serialize_model_(ostr);
  feasst_serialize_version(6798, ostr);
}

TwoBodyTable::TwoBodyTable(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(6798 == version, version);
}

}  // namespace feasst
