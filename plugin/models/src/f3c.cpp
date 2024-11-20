#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"
#include "configuration/include/physical_constants.h"
#include "models/include/f3c.h"

namespace feasst {

FEASST_MAPPER(F3C,);

F3C::F3C(argtype * args) {
  class_name_ = "F3C";
  Asc_ = dble("Asc", args, 1.);
}
F3C::F3C(argtype args) : F3C(&args) {
  feasst_check_all_used(args);
}

void F3C::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(2307, ostr);
  feasst_serialize(conversion_factor_, ostr);
  feasst_serialize(Asc_, ostr);
}

F3C::F3C(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2307, "unrecognized version: " << version);
  feasst_deserialize(&conversion_factor_, istr);
  feasst_deserialize(&Asc_, istr);
}

void F3C::precompute(const ModelParams& existing) {
  Model::precompute(existing);
  conversion_factor_ = existing.constants().charge_conversion();
}

double F3C::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double distance = std::sqrt(squared_distance);
  if (std::abs(distance) < NEAR_ZERO) {
    TRACE("near inf");
    return NEAR_INFINITY;
  }
  const double epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2],
    sigma = model_params.select(sigma_index()).mixed_values()[type1][type2],
    cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2],
    mixed_charge = model_params.select(charge_index()).mixed_values()[type1][type2],
    sig_inv_cut = sigma/cutoff,
    sig_inv_cut3 = sig_inv_cut*sig_inv_cut*sig_inv_cut,
    sig_inv_cut6 = sig_inv_cut3*sig_inv_cut3,
    const_shift = sig_inv_cut6*(Asc_*sig_inv_cut6 - 2.),
    linear_shift = -12.*sig_inv_cut6*(Asc_*sig_inv_cut6 - 1.)/cutoff,
    sig_inv_dis = sigma/distance,
    sig_inv_dis3 = sig_inv_dis*sig_inv_dis*sig_inv_dis,
    sig_inv_dis6 = sig_inv_dis3*sig_inv_dis3,
    en_lj = epsilon*(sig_inv_dis6*(Asc_*sig_inv_dis6 - 2.) - const_shift - linear_shift*(distance - cutoff)),
    const_shift_q = 1./cutoff,
    linear_shift_q = -const_shift_q*const_shift_q,
    en_q = mixed_charge*conversion_factor_*(1./distance - const_shift_q - linear_shift_q*(distance - cutoff));
  INFO("en_lj " << en_lj);
  INFO("en_q " << en_q);
  return en_lj + en_q;
}

}  // namespace feasst
