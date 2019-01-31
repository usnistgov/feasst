
#ifndef FEASST_EWALD_MODEL_CHARGE_SCREENED_H_
#define FEASST_EWALD_MODEL_CHARGE_SCREENED_H_

#include "core/include/model_two_body.h"
#include "core/include/constants.h"
#include "core/include/physical_constants.h"

namespace feasst {

/**
  Assumes units:
    1. length: Angstroms
    2. energy: kJ/mol
    3. charge: elementary
 */
class ModelChargeScreened : public ModelTwoBody {
 public:
  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    const double mixed_charge = model_params.mixed_charge()[type1][type2];
    const double distance = std::sqrt(squared_distance);
    const double alpha = model_params.property("alpha");
    return mixed_charge*charge_conversion*erfc(alpha*distance)/distance;
  }

  virtual ~ModelChargeScreened() {}
};

}  // namespace feasst

#endif  // FEASST_EWALD_MODEL_CHARGE_SCREENED_H_
