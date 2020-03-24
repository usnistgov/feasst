
#ifndef FEASST_SYSTEM_MODEL_TWO_BODY_H_
#define FEASST_SYSTEM_MODEL_TWO_BODY_H_

#include "system/include/model.h"
#include "system/include/visit_model.h"

namespace feasst {

class ModelTwoBody : public Model {
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
  double compute(
      const ModelParams& model_params,
      const Select& selection,
      const int group_index,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, model_params, selection, config, group_index);
    return visitor->energy();
  }
  double compute(
      const Select& selection,
      const int group_index,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, selection, config, group_index);
    return visitor->energy();
  }
  double compute(
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, model_params, selection, config, 0);
    return visitor->energy();
  }
  double compute(
      const Select& selection,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, selection, config, 0);
    return visitor->energy();
  }
  virtual ~ModelTwoBody() {}
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_TWO_BODY_H_
