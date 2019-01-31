
#ifndef FEASST_EXAMPLE_MODEL_EXAMPLE_H_
#define FEASST_EXAMPLE_MODEL_EXAMPLE_H_

#include "core/include/model_two_body.h"

namespace feasst {

/**
  Documentation for new class.
 */
class ModelExample : public ModelTwoBody {
 public:
  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    return 0;
  }
  virtual ~ModelExample() {}

 private:
};

}  // namespace feasst

#endif  // FEASST_EXAMPLE_MODEL_EXAMPLE_H_
