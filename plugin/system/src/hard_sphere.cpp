#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"
#include "system/include/hard_sphere.h"

namespace feasst {

FEASST_MAPPER(HardSphere,);

void HardSphere::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(607, ostr);
}

HardSphere::HardSphere(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(607 == version, version);
}

double HardSphere::energy(
  const double squared_distance,
  const int type1,
  const int type2,
  const ModelParams& model_params) {
  const double& sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  TRACE("sigma " << sigma);
  TRACE("sigma2 " << sigma*sigma);
  TRACE("r2 " << squared_distance);
  if (squared_distance <= sigma*sigma) {
    TRACE("near inf");
    return NEAR_INFINITY;
  }
  return 0.;
}

}  // namespace feasst
