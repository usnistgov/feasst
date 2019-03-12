
#ifndef FEASST_CORE_MODEL_LJ_SHIFT_H_
#define FEASST_CORE_MODEL_LJ_SHIFT_H_

#include "core/include/model_lj.h"

namespace feasst {

/**
  The Lennard-Jones potential, \f$U_{LJ}\f$ is described in ModelLJ.
  In this class, the \f$\alpha\f$ parameter may be varied (default: 6).
 */
class ModelLJAlpha : public ModelTwoBody {
 public:
  ModelLJAlpha() {
    set_alpha();
  }

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

  const double& hard_sphere_threshold_sq() const {
    return hard_sphere_threshold_sq_; }

  /// Return the derivative in the potential energy with respect to distance.
  double du_dr(
      const double distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const {
    const double epsilon = model_params.mixed_epsilon()[type1][type2];
    const double sigma = model_params.mixed_sigma()[type1][type2];
    const double rinv = sigma/distance;
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

  /// Initialize WCA cutoff
  void set_wca(const int site_type1, const int site_type2,
      ModelParams * params) {
    const double sigma = params->mixed_sigma()[site_type1][site_type2];
    const double r_wca = pow(2, 1./alpha())*sigma;
    params->set("cutoff", site_type1, site_type2, r_wca);
  }

  virtual ~ModelLJAlpha() {}

 private:
  double hard_sphere_threshold_sq_ = 0.2*0.2;
  double alpha_;
};

class EnergyAtCutoff : public ModelParam {
 public:
  void set_model(ModelLJAlpha * model) {
    model_ = model; }

  double compute(const int type1, const int type2,
      const ModelParams& model_params) override {
    const double cutoff = model_params.mixed_cutoff()[type1][type2];
    return model_->energy_without_shift(cutoff*cutoff, type1, type2, model_params);
  }

 private:
  ModelLJAlpha * model_;
};

/**
  The Lennard-Jones potential, \f$U_{LJ}\f$ is described in ModelLJ.
  This class implements the cut and shifted (CS) version which ensures \f$U(r_c)=0\f$.

  \f$ U_{LJ}^{CS}(r) = \left\{
    \begin{array}{lr}
      U_{LJ}(r) - U_{LJ}(r_c) & : r < r_c \\
      0 & r \ge r_c
    \end{array}
  \right. \f$
 */
class ModelLJCutShift : public ModelLJAlpha {
 public:
  // HWH some issues with this implementation include
  // - what if model params change, or is defined different by a special potential
  // - how to simplify user interface
  /// Precompute the shift factor for optimization, given existing model parameters.
  void precompute(const ModelParams& existing) override {
    precomputed_ = true;
    shift_.set_model(this); // note the model is used here for the computation
    shift_.set_param(existing);
    shift_.set_model(NULL); // remove model immediately
  }

  // Same as base class except with an error check.
  void set_alpha(const double alpha) override {
    ModelLJAlpha::set_alpha(alpha);
    ASSERT(!precomputed_, "shift depends on alpha. Set alpha first");
  }

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
    const double rinv_alpha = pow(rinv2, 0.5*alpha());
    const double en = 4.*epsilon*rinv_alpha*(rinv_alpha - 1.);
    const double shift = shift_.mixed_values()[type1][type2];
    return en - shift;
  }

  virtual ~ModelLJCutShift() {}

private:
  EnergyAtCutoff shift_;
  bool precomputed_ = false;
};

//  /**
//    The Weeks-Chandler-Andersen potential is the same as ModelLJCutShift, except
//    that the shifted cutoff occurs precisely at the potential minimum.
//
//    \f$ r_c = 2^{1/\alpha}\sigma \f$
//   */
//  class ModelWCA : public ModelLJCutShift {
//   public:
//    /// Precompute as in ModelLJCutShift, but also check that the cutoff is
//    /// consistent.
//    void precompute(const ModelParams& existing) override {
//      ModelLJCutShift::precompute(existing);
//
//      // check cutoff consistency
//      for (int type1 = 0; type1 < existing.size(); ++type1) {
//        for (int type2 = 0; type2 < existing.size(); ++type2) {
//          const double sigma = existing.mixed_sigma()[type1][type2];
//          const double cutoff = existing.mixed_cutoff()[type1][type2];
//          const double r_wca = pow(2, 1./alpha())*sigma;
//          ASSERT(std::abs(cutoff - r_wca) < NEAR_ZERO,
//            "the cutoff(" << cutoff << ") is not consistent with WCA");
//        }
//      }
//    }
//  };

class EnergyDerivAtCutoff : public ModelParam {
 public:
  void set_model(ModelLJAlpha * model) {
    model_ = model; }

  double compute(const int type1, const int type2,
      const ModelParams& model_params) override {
    const double cutoff = model_params.mixed_cutoff()[type1][type2];
    return model_->du_dr(cutoff, type1, type2, model_params);
  }

 private:
  ModelLJAlpha * model_;
};

/**
  The Lennard-Jones potential, \f$U_{LJ}\f$ is described in ModelLJ.
  This class implements the force shfited (FS) version which ensures both
  \f$U(r_c)=0\f$ and \f$\left.\frac{\partial U}{\partial r}\right|_{r=r_c}=0\f$.

  \f$ U_{LJ}^{FS}(r) = \left\{
    \begin{array}{lr}
      U_{LJ}(r) - U_{LJ}(r_c) - \left.\frac{\partial U}{\partial r}\right|_{r=r_c} (r - r_c) & : r < r_c \\
      0 & r \ge r_c
    \end{array}
  \right. \f$
 */
class ModelLJForceShift : public ModelLJAlpha {
 public:
  /// Precompute the shift factor for optimization, given existing model parameters.
  void precompute(const ModelParams& existing) override {
    precomputed_ = true;
    shift_.set_model(this); // note the model is used here for the computation
    shift_.set_param(existing);
    shift_.set_model(NULL); // remove model immediately

    force_shift_.set_model(this);
    force_shift_.set_param(existing);
    force_shift_.set_model(NULL);
  }

  void set_alpha(const double alpha) override {
    ModelLJAlpha::set_alpha(alpha);
    ASSERT(!precomputed_, "shift depends on alpha. Set alpha first");
  }

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
    const double rinv_alpha = pow(rinv2, 0.5*alpha());
    const double en = 4.*epsilon*rinv_alpha*(rinv_alpha - 1.);
    const double distance = std::sqrt(squared_distance);
    const double shift = shift_.mixed_values()[type1][type2];
    const double force_shift = force_shift_.mixed_values()[type1][type2];
    const double cutoff = model_params.mixed_cutoff()[type1][type2];
    return en - shift - (distance - cutoff)*force_shift;
  }

  virtual ~ModelLJForceShift() {}

 private:
  EnergyAtCutoff shift_;
  EnergyDerivAtCutoff force_shift_;
  bool precomputed_ = false;
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_LJ_SHIFT_H_
