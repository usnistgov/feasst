#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/configuration.h"
#include "models/include/lennard_jones_cut_shift.h"

namespace feasst {

FEASST_MAPPER(LennardJonesCutShift,);

LennardJonesCutShift::LennardJonesCutShift(argtype * args)
  : LennardJonesAlpha(args) {
  class_name_ = "LennardJonesCutShift";
}
LennardJonesCutShift::LennardJonesCutShift(argtype args) : LennardJonesCutShift(&args) {
  feasst_check_all_used(args);
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

void LennardJonesCutShift::precompute(const Configuration& config) {
  LennardJonesAlpha::precompute(config);
  const ModelParams& existing = config.model_params();
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
