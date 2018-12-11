
#ifndef FEASST_CORE_MODEL_HARD_SPHERE_H_
#define FEASST_CORE_MODEL_HARD_SPHERE_H_

#include "core/include/model_two_body.h"
#include "core/include/constants.h"

namespace feasst {

class ModelHardSphere : public ModelTwoBody {
 public:
  double evaluate(const Position &relative,
                  const Site& site1,
                  const Site& site2,
                  const ModelParams& model_params) const {
    const double squared_distance = relative.squared_distance();
    const int type1 = site1.type();
    const int type2 = site2.type();
    const double& sigma = (*model_params.mixed_sigma())[type1][type2];
    // std::cout << "sig " << sigma << " sqd " << squared_distance << std::endl;
    if (squared_distance <= sigma*sigma) {
      return NEAR_INFINITY;
    }
    return 0.;
  }
  virtual ~ModelHardSphere() {}
 private:
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_HARD_SPHERE_H_
