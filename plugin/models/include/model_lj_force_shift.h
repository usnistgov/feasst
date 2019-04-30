
#ifndef FEASST_MODELS_MODEL_LJ_FORCE_SHIFT_H_
#define FEASST_MODELS_MODEL_LJ_FORCE_SHIFT_H_

#include "models/include/model_lj_cut_shift.h"

namespace feasst {

class EnergyDerivAtCutoff : public ModelParam {
 public:
  EnergyDerivAtCutoff() : ModelParam() {}

  void set_model(ModelLJAlpha * model) {
    model_ = model; }

  double compute(const int type1, const int type2,
      const ModelParams& model_params) override {
    const double cutoff = model_params.mixed_cutoff()[type1][type2];
    return model_->du_dr(cutoff, type1, type2, model_params);
  }

  EnergyDerivAtCutoff(std::istream& istr) : ModelParam(istr) {}
 private:
  /// temporary pointer for use with compute without requiring new arguments.
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
  ModelLJForceShift() {}

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

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelLJForceShift>(istr); }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    serialize_model_lj_alpha_(ostr);
    feasst_serialize_version(923, ostr);
    shift_.serialize(ostr);
    force_shift_.serialize(ostr);
    feasst_serialize(precomputed_, ostr);
  }

  ModelLJForceShift(std::istream& istr) : ModelLJAlpha(istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(923 == version, version);
    shift_ = EnergyAtCutoff(istr);
    force_shift_ = EnergyDerivAtCutoff(istr);
    feasst_deserialize(&precomputed_, istr);
  }

  virtual ~ModelLJForceShift() {}

 private:
  const std::string class_name_ = "ModelLJForceShift";
  EnergyAtCutoff shift_;
  EnergyDerivAtCutoff force_shift_;
  bool precomputed_ = false;
};

inline std::shared_ptr<ModelLJForceShift> MakeModelLJForceShift() {
  return std::make_shared<ModelLJForceShift>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_MODEL_LJ_FORCE_SHIFT_H_
