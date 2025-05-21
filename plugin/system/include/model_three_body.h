
#ifndef FEASST_SYSTEM_MODEL_THREE_BODY_H_
#define FEASST_SYSTEM_MODEL_THREE_BODY_H_

#include "system/include/model.h"

namespace feasst {

class ModelTwoBody;
class VisitModel;

class ModelThreeBody : public Model {
 public:
  ModelThreeBody();

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
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) override;

  int num_body() const override { return 3; }

  ModelTwoBody * get_two_body() const;
  void set_two_body(std::unique_ptr<ModelTwoBody> two_body);

  virtual ~ModelThreeBody();
  explicit ModelThreeBody(std::istream& istr);

 private:
  std::unique_ptr<ModelTwoBody> two_body_;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_THREE_BODY_H_
