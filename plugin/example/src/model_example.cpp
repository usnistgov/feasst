
#include "example/include/model_example.h"

namespace feasst {

double ModelExample::evaluate(const Position &relative,
                              const Site& site1,
                              const Site& site2,
                              const ModelParams& model_params) const {
  const double squared_distance = relative.squared_distance();
  const int& type1 = site1.type();
  const int& type2 = site2.type();
  const double& cutoff = (*model_params.mixed_cutoff())[type1][type2];
  if (squared_distance <= cutoff*cutoff) {
    const double& sigma = (*model_params.mixed_sigma())[type1][type2];
    const double& epsilon = (*model_params.mixed_epsilon())[type1][type2];
    return 0*sigma*epsilon;
  }
  return 0.;
}

}  // namespace feasst
