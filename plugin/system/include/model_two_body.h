
#ifndef FEASST_SYSTEM_MODEL_TWO_BODY_H_
#define FEASST_SYSTEM_MODEL_TWO_BODY_H_

#include "system/include/model.h"

namespace feasst {

class VisitModel;

class ModelTwoBody : public Model {
 public:
  ModelTwoBody() {}

  double compute(
    const ModelParams& model_params,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) override;

  double compute(
    const int group_index,
    Configuration * config,
    VisitModel * visitor) override;

  double compute(
    const ModelParams& model_params,
    Configuration * config,
    VisitModel * visitor) override;

  double compute(
    Configuration * config,
    VisitModel * visitor) override;

  double compute(
    const ModelParams& model_params,
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) override;

  double compute(
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) override;

  double compute(
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) override;

  double compute(
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) override;

  int num_body() const override { return 2; }
  virtual ~ModelTwoBody() {}
  explicit ModelTwoBody(std::istream& istr) : Model(istr) {}
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_TWO_BODY_H_
