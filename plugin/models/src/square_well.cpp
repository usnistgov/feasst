#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"
#include "models/include/square_well.h"

namespace feasst {

FEASST_MAPPER(SquareWell,);

void SquareWell::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(553, ostr);
}

SquareWell::SquareWell(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 553, "unrecognized version: " << version);
}

double SquareWell::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double& sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  const double& epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2];
  if (squared_distance <= sigma*sigma) {
    TRACE("squared_distance " << squared_distance);
    return NEAR_INFINITY;
  }
  return -epsilon;
}

}  // namespace feasst
