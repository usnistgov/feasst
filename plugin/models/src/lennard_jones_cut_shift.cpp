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

LennardJonesCutShift::LennardJonesCutShift(argtype * args)
  : LennardJonesAlpha(args) {
  class_name_ = "LennardJonesCutShift";
}
LennardJonesCutShift::LennardJonesCutShift(argtype args) : LennardJonesCutShift(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void LennardJonesCutShift::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_lennard_jones_alpha_(ostr);
  feasst_serialize_version(644, ostr);
  shift_.serialize(ostr);
}

LennardJonesCutShift::LennardJonesCutShift(std::istream& istr)
  : LennardJonesAlpha(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(644 == version, version);
  shift_ = EnergyAtCutOff(istr);
}

void LennardJonesCutShift::precompute(const ModelParams& existing) {
  LennardJonesAlpha::precompute(existing);
  shift_.set_model(this); // note the model is used here for the computation
  shift_.set_param(existing);
  shift_.set_model(NULL); // remove model immediately
}

double LennardJonesCutShift::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double en = LennardJonesAlpha::energy(squared_distance, type1, type2, model_params);
  const double shift = shift_.mixed_values()[type1][type2];
  TRACE("en " << MAX_PRECISION << en);
  TRACE("shift " << MAX_PRECISION << shift);
  return en - shift;
}

}  // namespace feasst
