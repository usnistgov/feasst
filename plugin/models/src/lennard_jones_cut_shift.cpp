#include <cmath>
#include "models/include/lennard_jones_cut_shift.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"

namespace feasst {

class MapLennardJonesCutShift {
 public:
  MapLennardJonesCutShift() {
    LennardJonesCutShift().deserialize_map()["LennardJonesCutShift"] =
      MakeLennardJonesCutShift();
  }
};

static MapLennardJonesCutShift map_lennard_jones_cut_shift_ =
  MapLennardJonesCutShift();

LennardJonesCutShift::LennardJonesCutShift(argtype args)
  : LennardJonesAlpha(&args) {
  class_name_ = "LennardJonesCutShift";
}

void LennardJonesCutShift::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_lennard_jones_alpha_(ostr);
  feasst_serialize_version(644, ostr);
  shift_.serialize(ostr);
  ostr << precomputed_ << " ";
}

LennardJonesCutShift::LennardJonesCutShift(std::istream& istr)
  : LennardJonesAlpha(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(644 == version, version);
  shift_ = EnergyAtCutoff(istr);
  istr >> precomputed_;
}

void LennardJonesCutShift::precompute(const ModelParams& existing) {
  precomputed_ = true;
  shift_.set_model(this); // note the model is used here for the computation
  shift_.set_param(existing);
  shift_.set_model(NULL); // remove model immediately
}

void LennardJonesCutShift::set_wca(const int site_type1, const int site_type2,
    ModelParams * params) {
  const double sigma = params->mixed_sigma()[site_type1][site_type2];
  const double r_wca = std::pow(2, 1./alpha())*sigma;
  params->set("cutoff", site_type1, site_type2, r_wca);
}

double LennardJonesCutShift::energy(
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
  const double shift = shift_.mixed_values()[type1][type2];
  return en - shift;
}

}  // namespace feasst
