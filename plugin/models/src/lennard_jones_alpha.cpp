#include <cmath>
#include "models/include/lennard_jones_alpha.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"

namespace feasst {

LennardJonesAlpha::LennardJonesAlpha(const argtype& args) : LennardJones(args) {
  alpha_ = args_.key("alpha").dflt("6").dble();
  class_name_ = "LennardJonesAlpha";
}

class MapLennardJonesAlpha {
 public:
  MapLennardJonesAlpha() {
    LennardJonesAlpha().deserialize_map()["LennardJonesAlpha"] = MakeLennardJonesAlpha();
  }
};

static MapLennardJonesAlpha map_lennard_jones_alpha_ = MapLennardJonesAlpha();

void LennardJonesAlpha::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_lennard_jones_alpha_(ostr);
}

void LennardJonesAlpha::serialize_lennard_jones_alpha_(
    std::ostream& ostr) const {
  serialize_lennard_jones_(ostr);
  feasst_serialize_version(713, ostr);
  feasst_serialize(alpha_, ostr);
}

LennardJonesAlpha::LennardJonesAlpha(std::istream& istr) : LennardJones(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(713 == version, version);
  feasst_deserialize(&alpha_, istr);
}

double LennardJonesAlpha::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double sigma = model_params.mixed_sigma()[type1][type2];
  const double sigma_squared = sigma*sigma;
  if (squared_distance < hard_sphere_threshold_sq()*sigma_squared) {
    return NEAR_INFINITY;
  }
  const double epsilon = model_params.mixed_epsilon()[type1][type2];
  const double rinv2 = sigma_squared/squared_distance;
  const double rinv_alpha = std::pow(rinv2, 0.5*alpha_);
  return 4.*epsilon*rinv_alpha*(rinv_alpha - 1.);
}

double LennardJonesAlpha::du_dr(
    const double distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) const {
  const double epsilon = model_params.mixed_epsilon()[type1][type2];
  const double sigma = model_params.mixed_sigma()[type1][type2];
  const double rinv = sigma/distance;
  if (sigma == 0) {
    return 0.;
  }
  return 4.*epsilon/sigma*alpha()*(-2*std::pow(rinv, 2*alpha() + 1)
                                    + std::pow(rinv, alpha() + 1));
}

}  // namespace feasst
