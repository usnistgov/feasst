#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/table.h"
#include "configuration/include/physical_constants.h"
#include "configuration/include/model_params.h"
#include "configuration/include/physical_constants.h"
#include "charge/include/charge_screened.h"

namespace feasst {

class MapChargeScreened {
 public:
  MapChargeScreened() {
    ChargeScreened().deserialize_map()["ChargeScreened"] = MakeChargeScreened();
  }
};

static MapChargeScreened map_charge_screened_ = MapChargeScreened();

ChargeScreened::ChargeScreened(argtype * args) {
  class_name_ = "ChargeScreened";
  erfc_table_size_ = integer("erfc_table_size", args, 0);
  const double hs_thres = dble("hard_sphere_threshold", args, 0.2);
  hard_sphere_threshold_sq_ = hs_thres*hs_thres;
}
ChargeScreened::ChargeScreened(argtype args) : ChargeScreened(&args) {
  feasst_check_all_used(args);
}

void ChargeScreened::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(4616, ostr);
  feasst_serialize(alpha_, ostr);
  feasst_serialize(conversion_factor_, ostr);
  feasst_serialize(hard_sphere_threshold_sq_, ostr);
  feasst_serialize(erfc_table_size_, ostr);
  feasst_serialize(erfc_, ostr);
}

ChargeScreened::ChargeScreened(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4616, "unrecognized verison: " << version);
  feasst_deserialize(&alpha_, istr);
  feasst_deserialize(&conversion_factor_, istr);
  feasst_deserialize(&hard_sphere_threshold_sq_, istr);
  feasst_deserialize(&erfc_table_size_, istr);
  // HWH for unknown reasons, this function template does not work.
  // feasst_deserialize_fstdr(erfc_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      erfc_ = std::make_shared<Table1D>(istr);//erfc_->deserialize(istr);
    }
  }
}

double ChargeScreened::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double mixed_charge = model_params.select(charge_index()).mixed_values()[type1][type2];
  if (squared_distance < hard_sphere_threshold_sq_) {
    return NEAR_INFINITY;
  } else if (erfc_) {
    const double mixed_max_cutoff = model_params.select(cutoff_index()).mixed_max();
    const double z = squared_distance/mixed_max_cutoff/mixed_max_cutoff;
    const double erffac = erfc_->forward_difference_interpolation(z);
    //const double erffac = erfc_->linear_interpolation(z);
    TRACE("erffac " << erffac);
    TRACE("U " << mixed_charge*conversion_factor_*erffac);
    return mixed_charge*conversion_factor_*erffac;
  } else {
    const double distance = std::sqrt(squared_distance);
    const double en = mixed_charge*conversion_factor_*std::erfc(alpha_*distance)/distance;
    //TRACE("mixed_charge " << mixed_charge);
    //TRACE("conversion_factor_ " << conversion_factor_);
    TRACE("en " << en);
    return en;
  }
}

void ChargeScreened::precompute(const ModelParams& existing) {
  Model::precompute(existing);
  alpha_ = existing.property("alpha");
  conversion_factor_ = existing.constants().charge_conversion();
  init_erfc_(existing.select("cutoff").mixed_max());
}

void ChargeScreened::init_erfc_(const double cutoff) {
  if (erfc_table_size_ > 0) {
    erfc_ = MakeTable1D({{"num", str(erfc_table_size_)}});
    for (int bin = 0; bin < erfc_->num(); ++bin) {
      const double z = erfc_->bin_to_value(bin);
      const double x = std::sqrt(z)*cutoff;
      erfc_->set_data(bin, std::erfc(alpha_*x)/x);
      if (bin == 1) {
        ASSERT(x*x < hard_sphere_threshold_sq_,
          "erfc_table_size: " << erfc_table_size_ << " leads to a spacing " <<
          "of the first bin as " << x << ". This is larger than the " <<
          "hard_sphere_threshold: " << std::sqrt(hard_sphere_threshold_sq_) <<
          ". Either increase this threadhold or the table size.");
      }
    }
  }
}

}  // namespace feasst
