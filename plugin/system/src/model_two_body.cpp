#include "system/include/visit_model.h"
#include "system/include/model_two_body.h"

namespace feasst {

ModelTwoBody::ModelTwoBody() {}
ModelTwoBody::~ModelTwoBody() {}

double ModelTwoBody::compute(
    const ModelParams& model_params,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, config, group_index);
  return visitor->energy();
}

double ModelTwoBody::compute(
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, config, group_index);
  return visitor->energy();
}

double ModelTwoBody::compute(
    const ModelParams& model_params,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, config, 0);
  return visitor->energy();
}

double ModelTwoBody::compute(
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, config, 0);
  return visitor->energy();
}

double ModelTwoBody::compute(
    const ModelParams& model_params,
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, selection, config, group_index);
  return visitor->energy();
}

double ModelTwoBody::compute(
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, selection, config, group_index);
  return visitor->energy();
}

double ModelTwoBody::compute(
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, selection, config, 0);
  return visitor->energy();
}

double ModelTwoBody::compute(
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, selection, config, 0);
  return visitor->energy();
}

}  // namespace feasst
