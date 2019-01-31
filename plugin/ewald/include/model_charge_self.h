
#ifndef FEASST_EWALD_MODEL_CHARGE_SELF_H_
#define FEASST_EWALD_MODEL_CHARGE_SELF_H_

#include "core/include/model_one_body.h"
#include "core/include/constants.h"
#include "core/include/physical_constants.h"

namespace feasst {

/**
  Assumes units:
    1. length: Angstroms
    2. energy: kJ/mol
    3. charge: elementary
 */
class ModelChargeSelf : public ModelOneBody {
 public:
  double energy(
      const Site& site,
      const Configuration * config,
      const ModelParams& model_params) const {
    const int type = site.type();
    const double charge = model_params.charge().value(type);
    const double alpha = model_params.property("alpha");
    return -charge*charge*charge_conversion*alpha/std::sqrt(PI);
  }

  virtual ~ModelChargeSelf() {}
};

}  // namespace feasst

#endif  // FEASST_EWALD_MODEL_CHARGE_SELF_H_
