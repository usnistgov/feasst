
#ifndef FEASST_CORE_MODEL_HARD_SPHERE_H_
#define FEASST_CORE_MODEL_HARD_SPHERE_H_

#include "core/include/model_two_body.h"
#include "core/include/constants.h"

namespace feasst {

class ModelHardSphere : public ModelTwoBody {
 public:
  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) const override {
    const double& sigma = model_params.mixed_sigma()[type1][type2];
    if (squared_distance <= sigma*sigma) {
      return NEAR_INFINITY;
    }
    return 0.;
  }
  virtual ~ModelHardSphere() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_HARD_SPHERE_H_
