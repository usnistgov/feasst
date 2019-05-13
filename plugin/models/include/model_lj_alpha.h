
#ifndef FEASST_MODELS_MODEL_LJ_ALPHA_H_
#define FEASST_MODELS_MODEL_LJ_ALPHA_H_

#include "system/include/model_lj.h"

namespace feasst {

/**
  The Lennard-Jones potential, \f$U_{LJ}\f$ is described in ModelLJ.
  In this class, the \f$\alpha\f$ parameter may be varied (default: 6).
 */
class ModelLJAlpha : public ModelLJ {
 public:
  ModelLJAlpha();

  /// Set the value of \f$\alpha\f$. The default value is 6.
  virtual void set_alpha(const double alpha = 6) { alpha_ = alpha; }

  /// Return the value of alpha.
  const double& alpha() const { return alpha_; }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    const double sigma = model_params.mixed_sigma()[type1][type2];
    const double sigma_squared = sigma*sigma;
    if (squared_distance < hard_sphere_threshold_sq()*sigma_squared) {
      return NEAR_INFINITY;
    }
    const double epsilon = model_params.mixed_epsilon()[type1][type2];
    const double rinv2 = sigma_squared/squared_distance;
    const double rinv_alpha = pow(rinv2, 0.5*alpha_);
    return 4.*epsilon*rinv_alpha*(rinv_alpha - 1.);
  }

  /// Return the derivative in the potential energy with respect to distance.
  double du_dr(
      const double distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const {
    const double epsilon = model_params.mixed_epsilon()[type1][type2];
    const double sigma = model_params.mixed_sigma()[type1][type2];
    const double rinv = sigma/distance;
    if (sigma == 0) {
      return 0.;
    }
    return 4.*epsilon/sigma*alpha()*(-2*pow(rinv, 2*alpha() + 1)
                                      + pow(rinv, alpha() + 1));
  }

  // Derived classes may call the ModelLJAlpha energy directly.
  double energy_without_shift(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const {
    return ModelLJAlpha::energy(squared_distance, type1, type2, model_params);
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelLJAlpha>(istr);
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    serialize_model_lj_alpha_(ostr);
  }

  ModelLJAlpha(std::istream& istr) : ModelLJ(istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(713 == version, version);
    feasst_deserialize(&alpha_, istr);
  }

  virtual ~ModelLJAlpha() {}

 protected:
  void serialize_model_lj_alpha_(std::ostream& ostr) const {
    serialize_model_lj_(ostr);
    feasst_serialize_version(713, ostr);
    feasst_serialize(alpha_, ostr);
  }

 private:
  const std::string class_name_ = "ModelLJAlpha";
  double alpha_;
};

inline std::shared_ptr<ModelLJAlpha> MakeModelLJAlpha() {
  return std::make_shared<ModelLJAlpha>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_MODEL_LJ_ALPHA_H_
