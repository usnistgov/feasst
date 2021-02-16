#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"
#include "system/include/lennard_jones.h"

namespace feasst {

class MapLennardJones {
 public:
  MapLennardJones() {
    LennardJones().deserialize_map()["LennardJones"] = MakeLennardJones();
  }
};

static MapLennardJones mapper_ = MapLennardJones();

LennardJones::LennardJones(argtype args) : LennardJones(&args) {
  check_all_used(args); }
LennardJones::LennardJones(argtype * args) {
  class_name_ = "LennardJones";
  const double thres = dble("hard_sphere_threshold", args, 0.2);
  hard_sphere_threshold_sq_ = thres*thres;
}

void LennardJones::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_lennard_jones_(ostr);
}

void LennardJones::serialize_lennard_jones_(std::ostream& ostr) const {
  feasst_serialize_version(763, ostr);
  feasst_serialize(hard_sphere_threshold_sq_, ostr);
}

LennardJones::LennardJones(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(763 == version, version);
  feasst_deserialize(&hard_sphere_threshold_sq_, istr);
}

double LennardJones::hard_sphere_threshold() const {
  return std::sqrt(hard_sphere_threshold_sq_);
}

double LennardJones::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  TRACE("squared_distance " << squared_distance);
  TRACE("type1 " << type1);
  TRACE("type2 " << type2);
  const double sigma = model_params.mixed_sigma()[type1][type2];
  const double sigma_squared = sigma*sigma;
  if (squared_distance < hard_sphere_threshold_sq_*sigma_squared) {
    return NEAR_INFINITY;
  }
  const double epsilon = model_params.mixed_epsilon()[type1][type2];
  const double rinv2 = sigma_squared/squared_distance;
  const double rinv6 = rinv2*rinv2*rinv2;
  const double en = 4.*epsilon*rinv6*(rinv6 - 1.);
  TRACE("en " << en);
  return en;
}

}  // namespace feasst
