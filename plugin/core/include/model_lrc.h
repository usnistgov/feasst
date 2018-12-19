
#ifndef FEASST_CORE_MODEL_LRC_H_
#define FEASST_CORE_MODEL_LRC_H_

#include "core/include/model_one_body.h"

namespace feasst {

class ModelLRC : public ModelOneBody {
 public:
  double evaluate(const Site& site,
                  const Configuration& config,
                  const ModelParams& model_params) const {
    const int type = site.type();
    const double epsilon = model_params.epsilon().value(type);
    const double sigma = model_params.sigma().value(type);
    const double cut_distance = model_params.cutoff().value(type);
    const double prefactor = epsilon*(8./3.)*PI*pow(sigma, 3)*((1./3.)*pow(sigma/cut_distance, 9) - pow(sigma/cut_distance, 3));
    return double(config.num_particles()/config.domain().volume())*prefactor;
  }

  virtual ~ModelLRC() {}
 private:
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_LRC_H_
