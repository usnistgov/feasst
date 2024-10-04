#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "configuration/include/model_params.h"
#include "system/include/cutoff_outer.h"

namespace feasst {

FEASST_MAPPER(CutoffOuter,);

void CutoffOuter::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(2498, ostr);
}

CutoffOuter::CutoffOuter(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2498, "mismatch version: " << version);
}

}  // namespace feasst
