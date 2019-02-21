
#ifndef FEASST_CORE_MODEL_SQUARE_WELL_H_
#define FEASST_CORE_MODEL_SQUARE_WELL_H_

#include "core/include/model_two_body.h"
#include "core/include/constants.h"

namespace feasst {

class ModelSquareWell : public ModelTwoBody {
 public:
  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) const override {
    const double& sigma = model_params.mixed_sigma()[type1][type2];
    const double& epsilon = model_params.mixed_epsilon()[type1][type2];
    if (squared_distance <= sigma*sigma) {
      return NEAR_INFINITY;
    }
    return -epsilon;
  }
  virtual ~ModelSquareWell() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_SQUARE_WELL_H_
