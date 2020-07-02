#include "configuration/include/site.h"
#include "system/include/model_one_body.h"

namespace feasst {

double ModelOneBody::energy_no_wrap(
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const {
  return energy(site.position(), site, config, model_params);
}

}  // namespace feasst
