
#ifndef FEASST_CORE_MODEL_YUKAWA_H_
#define FEASST_CORE_MODEL_YUKAWA_H_

#include "core/include/model_two_body.h"
#include "core/include/constants.h"

namespace feasst {

/**
  The screened Coulomb (Yukawa) potential is given by
  \f$ U = \epsilon \exp(-\kappa (r/\sigma - 1)) / (r/\sigma) \f$
  Note that \f$\kappa\f$ is dimensionless in this definition.
 */
class ModelYukawa : public ModelTwoBody {
 public:
  ModelYukawa() { set_kappa(); }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    const double epsilon = model_params.mixed_epsilon()[type1][type2];
    const double sigma = model_params.mixed_sigma()[type1][type2];
    const double distance = sqrt(squared_distance);
    TRACE("epsilon " << epsilon << " distance " << distance << " kappa " << kappa_);
    return epsilon*exp(-kappa_*(distance/sigma - 1.))/(distance/sigma);
  }

  /// Set the value of the kappa parameter.
  void set_kappa(const double kappa = 1.) { kappa_ = kappa; }

  virtual ~ModelYukawa() {}

 private:
  double kappa_;
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_YUKAWA_H_
