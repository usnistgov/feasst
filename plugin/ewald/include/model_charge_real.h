
#ifndef FEASST_EWALD_MODEL_CHARGE_REAL_H_
#define FEASST_EWALD_MODEL_CHARGE_REAL_H_

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
class ModelChargeReal : public ModelTwoBody {
 public:
  /// HWH error check that alpha is set.
  void set_alpha(const double alpha) { alpha_ = alpha; }
  double alpha() const { return alpha_; }

  double evaluate(const Position &relative,
                  const Site& site1,
                  const Site& site2,
                  const ModelParams& model_params) const {
    const double squared_distance = relative.squared_distance();
    const int& type1 = site1.type();
    const int& type2 = site2.type();
    const double& cutoff = (*model_params.mixed_cutoff())[type1][type2];
    if (squared_distance <= cutoff*cutoff) {
      const double& mixed_charge = (*model_params.mixed_charge())[type1][type2];
      const double distance = std::sqrt(squared_distance);
      return mixed_charge*charge_conversion_*erfc(alpha_*distance)/distance;
    }
    return 0.;
  }

  virtual ~ModelChargeReal() {}

 private:
  double alpha_;
  const double charge_conversion_ = pow(elementary_charge, 2)/
              (4*PI*permitivity_vacuum*1e3/1e10/avogadro_constant);
};

}  // namespace feasst

#endif  // FEASST_EWALD_MODEL_CHARGE_REAL_H_
