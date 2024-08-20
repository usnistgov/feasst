#include "configuration/include/site.h"
#include "system/include/visit_model.h"
#include "system/include/model_one_body.h"

namespace feasst {

double ModelOneBody::energy_no_wrap(
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) {
  return energy(site.position(), site, config, model_params);
}

double ModelOneBody::compute(
    const ModelParams& model_params,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, config, group_index);
  return visitor->energy();
}

double ModelOneBody::compute(
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, config, group_index);
  return visitor->energy();
}

double ModelOneBody::compute(
    const ModelParams& model_params,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, config, 0);
  return visitor->energy();
}

double ModelOneBody::compute(
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, config, 0);
  return visitor->energy();
}

double ModelOneBody::compute(
    const ModelParams& model_params,
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, selection, config, group_index);
  return visitor->energy();
}

double ModelOneBody::compute(
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, selection, config, group_index);
  return visitor->energy();
}

double ModelOneBody::compute(
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, model_params, selection, config, 0);
  return visitor->energy();
}

double ModelOneBody::compute(
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) {
  visitor->compute(this, selection, config, 0);
  return visitor->energy();
}

}  // namespace feasst
