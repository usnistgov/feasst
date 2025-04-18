#include "utils/include/debug.h"
#include "system/include/visit_model.h"
#include "system/include/model_two_body.h"
#include "system/include/model_three_body.h"

namespace feasst {

double ModelThreeBody::compute(
    const ModelParams& model_params,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, config, group_index);
  return visitor->energy();
}

double ModelThreeBody::compute(
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, config, group_index);
  return visitor->energy();
}

double ModelThreeBody::compute(
    const ModelParams& model_params,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, config, 0);
  return visitor->energy();
}

double ModelThreeBody::compute(
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, config, 0);
  return visitor->energy();
}

double ModelThreeBody::compute(
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, selection, config, 0);
  return visitor->energy();
}

double ModelThreeBody::compute(
    const ModelParams& model_params,
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, selection, config, group_index);
  return visitor->energy();
}

double ModelThreeBody::compute(
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, selection, config, group_index);
  return visitor->energy();
}

double ModelThreeBody::compute(
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, selection, config, 0);
  return visitor->energy();
}

ModelTwoBody * ModelThreeBody::get_two_body() const {
  return two_body_.get();
}

void ModelThreeBody::set_two_body(std::unique_ptr<ModelTwoBody> two_body) {
  two_body_ = std::move(two_body);
}

}  // namespace feasst
