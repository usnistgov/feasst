
#ifndef FEASST_SYSTEM_MODEL_THREE_BODY_H_
#define FEASST_SYSTEM_MODEL_THREE_BODY_H_

#include "system/include/model.h"
#include "system/include/visit_model.h"
#include "configuration/include/particle.h"
#include "math/include/constants.h"

namespace feasst {

class ModelThreeBody : public Model {
 public:
  double compute(
      const ModelParams& model_params,
      const int group_index,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, model_params, config, group_index);
    return visitor->energy();
  }
  double compute(
      const int group_index,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, config, group_index);
    return visitor->energy();
  }
  double compute(
      const ModelParams& model_params,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, model_params, config, 0);
    return visitor->energy();
  }
  double compute(
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, config, 0);
    return visitor->energy();
  }
  virtual double energy(
      const int type1,
      const int type2,
      const int type3,
      const ModelParams& model_params) const = 0;
  virtual ~ModelThreeBody() {}
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_THREE_BODY_H_
