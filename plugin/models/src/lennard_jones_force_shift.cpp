#include <cmath>
#include "models/include/lennard_jones_force_shift.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"

namespace feasst {

class MapLennardJonesForceShift {
 public:
  MapLennardJonesForceShift() {
    LennardJonesForceShift().deserialize_map()["LennardJonesForceShift"] = MakeLennardJonesForceShift();
  }
};

static MapLennardJonesForceShift map_lennard_jones_force_shift_ = MapLennardJonesForceShift();

LennardJonesForceShift::LennardJonesForceShift(argtype args)
  : LennardJonesAlpha(&args) {
  class_name_ = "LennardJonesForceShift";
}

void LennardJonesForceShift::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_lennard_jones_alpha_(ostr);
  feasst_serialize_version(923, ostr);
  shift_.serialize(ostr);
  force_shift_.serialize(ostr);
  feasst_serialize(precomputed_, ostr);
}

LennardJonesForceShift::LennardJonesForceShift(std::istream& istr)
  : LennardJonesAlpha(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(923 == version, version);
  shift_ = EnergyAtCutoff(istr);
  force_shift_ = EnergyDerivAtCutoff(istr);
  feasst_deserialize(&precomputed_, istr);
}

void LennardJonesForceShift::precompute(const ModelParams& existing) {
  precomputed_ = true;
  shift_.set_model(this); // note the model is used here for the computation
  shift_.set_param(existing);
  shift_.set_model(NULL); // remove model immediately

  force_shift_.set_model(this);
  force_shift_.set_param(existing);
  force_shift_.set_model(NULL);
}

double LennardJonesForceShift::energy(
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
  const double rinv_alpha = std::pow(rinv2, 0.5*alpha());
  const double en = 4.*epsilon*rinv_alpha*(rinv_alpha - 1.);
  const double distance = std::sqrt(squared_distance);
  const double shift = shift_.mixed_values()[type1][type2];
  const double force_shift = force_shift_.mixed_values()[type1][type2];
  const double cutoff = model_params.mixed_cutoff()[type1][type2];
  return en - shift - (distance - cutoff)*force_shift;
}

}  // namespace feasst
