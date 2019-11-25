
#ifndef FEASST_MODELS_LENNARD_JONES_CUT_SHIFT_H_
#define FEASST_MODELS_LENNARD_JONES_CUT_SHIFT_H_

#include "models/include/lennard_jones_alpha.h"

namespace feasst {

/**
  The Lennard-Jones potential, \f$U_{LJ}\f$ is described in LennardJones.
  This class implements the cut and shifted (CS) version which ensures \f$U(r_c)=0\f$.

  \f$ U_{LJ}^{CS}(r) = \left\{
    \begin{array}{lr}
      U_{LJ}(r) - U_{LJ}(r_c) & : r < r_c \\
      0 & r \ge r_c
    \end{array}
  \right. \f$

 For a Weeks-Chandler-Anderson (WCA) potential, use this class.
 The cutoffs are computed by LennardJonesCutShift::set_wca.
 Thus, one workflow is as follows:

   ModelParams wca_params = configuration.model_params(); // copy model params
   LennardJonesCutShift::set_wca(type1, type2, wca_params);  // set cutoff
   wca->precompute(wca_params);   // compute shifts, etc
   Potential::set_model_params(wca_params); // use wca_params.

 */
class LennardJonesCutShift : public LennardJonesAlpha {
 public:
  LennardJonesCutShift() {}

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
    LennardJonesAlpha::set_alpha(alpha);
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
    return std::make_shared<LennardJonesCutShift>(istr); }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    serialize_lennard_jones_alpha_(ostr);
    feasst_serialize_version(644, ostr);
    shift_.serialize(ostr);
    ostr << precomputed_ << " ";
  }

  LennardJonesCutShift(std::istream& istr) : LennardJonesAlpha(istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(644 == version, version);
    shift_ = EnergyAtCutoff(istr);
    istr >> precomputed_;
  }

  virtual ~LennardJonesCutShift() {}

private:
  const std::string class_name_ = "LennardJonesCutShift";
  EnergyAtCutoff shift_;
  bool precomputed_ = false;
};

inline std::shared_ptr<LennardJonesCutShift> MakeLennardJonesCutShift() {
  return std::make_shared<LennardJonesCutShift>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_LENNARD_JONES_CUT_SHIFT_H_
