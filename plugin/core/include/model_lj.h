
#ifndef FEASST_CORE_MODEL_LJ_H_
#define FEASST_CORE_MODEL_LJ_H_

#include "core/include/model_two_body.h"
#include "core/include/constants.h"

namespace feasst {

class ModelLJ : public ModelTwoBody {
 public:
  double evaluate(const Position &relative,
                  const Site& site1,
                  const Site& site2,
                  const ModelParams& model_params) const {
    const double squared_distance = relative.squared_distance();
    const int type1 = site1.type();
    const int type2 = site2.type();
    const double cutoff = (*model_params.mixed_cutoff())[type1][type2];
    TRACE("squared dist: " << squared_distance);
    if (squared_distance <= cutoff*cutoff) {
      const double sigma = (*model_params.mixed_sigma())[type1][type2];
      const double sigma_squared = sigma*sigma;
      TRACE("sigsq: " << sigma_squared);
      if (squared_distance < hard_sphere_threshold_sq_*sigma_squared) {
        return NEAR_INFINITY;
      }
      const double epsilon = (*model_params.mixed_epsilon())[type1][type2];
      const double rinv2 = sigma_squared/squared_distance;
      const double rinv6 = rinv2*rinv2*rinv2;
      return 4.*epsilon*rinv6*(rinv6 - 1.);
    }
    return 0.;
  }

  virtual ~ModelLJ() {}

 private:
  double hard_sphere_threshold_sq_ = 0.2*0.2;
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_LJ_H_
