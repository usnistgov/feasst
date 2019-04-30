
#ifndef FEASST_MODELS_MODEL_LJ_CUT_SHIFT_H_
#define FEASST_MODELS_MODEL_LJ_CUT_SHIFT_H_

#include "models/include/model_lj_alpha.h"

namespace feasst {

class EnergyAtCutoff : public ModelParam {
 public:
  EnergyAtCutoff() : ModelParam() {}

  void set_model(ModelLJAlpha * model) {
    model_ = model; }

  double compute(const int type1, const int type2,
      const ModelParams& model_params) override {
    const double cutoff = model_params.mixed_cutoff()[type1][type2];
    return model_->energy_without_shift(cutoff*cutoff, type1, type2, model_params);
  }

  EnergyAtCutoff(std::istream& istr) : ModelParam(istr) {}

private:
  /// temporary pointer for use with compute without requiring new arguments.
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

 For a Weeks-Chandler-Anderson (WCA) potential, use this class.
 The cutoffs are computed by ModelLJCutShift::set_wca.
 Thus, one workflow is as follows:

   ModelParams wca_params = configuration.model_params(); // copy model params
   ModelLJCutShift::set_wca(type1, type2, wca_params);  // set cutoff
   wca->precompute(wca_params);   // compute shifts, etc
   Potential::set_model_params(wca_params); // use wca_params.

 */
class ModelLJCutShift : public ModelLJAlpha {
 public:
  ModelLJCutShift() {}

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

  /// Initialize WCA cutoff distances for types 1 and 2 in model parameters.
  void set_wca(const int site_type1, const int site_type2,
      ModelParams * params) {
    const double sigma = params->mixed_sigma()[site_type1][site_type2];
    const double r_wca = pow(2, 1./alpha())*sigma;
    params->set("cutoff", site_type1, site_type2, r_wca);
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelLJCutShift>(istr); }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    serialize_model_lj_alpha_(ostr);
    feasst_serialize_version(644, ostr);
    shift_.serialize(ostr);
    ostr << precomputed_ << " ";
  }

  ModelLJCutShift(std::istream& istr) : ModelLJAlpha(istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(644 == version, version);
    shift_ = EnergyAtCutoff(istr);
    istr >> precomputed_;
  }

  virtual ~ModelLJCutShift() {}

private:
  const std::string class_name_ = "ModelLJCutShift";
  EnergyAtCutoff shift_;
  bool precomputed_ = false;
};

inline std::shared_ptr<ModelLJCutShift> MakeModelLJCutShift() {
  return std::make_shared<ModelLJCutShift>();
}

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

}  // namespace feasst

#endif  // FEASST_MODELS_MODEL_LJ_CUT_SHIFT_H_
