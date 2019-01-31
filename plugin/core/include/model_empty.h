
#ifndef FEASST_CORE_MODEL_EMPTY_H_
#define FEASST_CORE_MODEL_EMPTY_H_

#include "core/include/model_one_body.h"

namespace feasst {

class ModelEmpty : public ModelOneBody {
 public:
  double energy(
      const Site& site,
      const Configuration * config,
      const ModelParams& model_params) const { return 0.; }
  virtual ~ModelEmpty() {}
 private:
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_EMPTY_H_
