#include "system/include/visit_model.h"
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

}  // namespace feasst
