
#ifndef FEASST_CORE_MODEL_LJ_H_
#define FEASST_CORE_MODEL_LJ_H_

#include "core/include/model_two_body.h"
#include "core/include/constants.h"

namespace feasst {

class ModelLJ : public ModelTwoBody {
 public:
  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    const double sigma = model_params.mixed_sigma()[type1][type2];
    const double sigma_squared = sigma*sigma;
    if (squared_distance < hard_sphere_threshold_sq_*sigma_squared) {
      return NEAR_INFINITY;
    }
    const double epsilon = model_params.mixed_epsilon()[type1][type2];
    const double rinv2 = sigma_squared/squared_distance;
    const double rinv6 = rinv2*rinv2*rinv2;
    return 4.*epsilon*rinv6*(rinv6 - 1.);
  }

  virtual ~ModelLJ() {}

 private:
  double hard_sphere_threshold_sq_ = 0.2*0.2;
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_LJ_H_
