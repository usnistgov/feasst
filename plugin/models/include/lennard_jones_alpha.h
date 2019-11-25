
#ifndef FEASST_MODELS_LENNARD_JONES_ALPHA_H_
#define FEASST_MODELS_LENNARD_JONES_ALPHA_H_

#include "system/include/lennard_jones.h"

namespace feasst {

/**
  The Lennard-Jones potential, \f$U_{LJ}\f$ is described in LennardJones.
  In this class, the \f$\alpha\f$ parameter may be varied (default: 6).
 */
class LennardJonesAlpha : public LennardJones {
 public:
  LennardJonesAlpha();

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

  // Derived classes may call the LennardJonesAlpha energy directly.
  double energy_without_shift(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const {
    return LennardJonesAlpha::energy(squared_distance, type1, type2, model_params);
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<LennardJonesAlpha>(istr);
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    serialize_lennard_jones_alpha_(ostr);
  }

  LennardJonesAlpha(std::istream& istr) : LennardJones(istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(713 == version, version);
    feasst_deserialize(&alpha_, istr);
  }

  virtual ~LennardJonesAlpha() {}

 protected:
  void serialize_lennard_jones_alpha_(std::ostream& ostr) const {
    serialize_lennard_jones_(ostr);
    feasst_serialize_version(713, ostr);
    feasst_serialize(alpha_, ostr);
  }

 private:
  const std::string class_name_ = "LennardJonesAlpha";
  double alpha_;
};

inline std::shared_ptr<LennardJonesAlpha> MakeLennardJonesAlpha() {
  return std::make_shared<LennardJonesAlpha>();
}

class EnergyAtCutoff : public ModelParam {
 public:
  EnergyAtCutoff() : ModelParam() {}

  void set_model(LennardJonesAlpha * model) {
    model_ = model; }

  double compute(const int type1, const int type2,
      const ModelParams& model_params) override {
    const double cutoff = model_params.mixed_cutoff()[type1][type2];
    return model_->energy_without_shift(cutoff*cutoff, type1, type2, model_params);
  }

  EnergyAtCutoff(std::istream& istr) : ModelParam(istr) {}

private:
  /// temporary pointer for use with compute without requiring new arguments.
  LennardJonesAlpha * model_;
};

class EnergyDerivAtCutoff : public ModelParam {
 public:
  EnergyDerivAtCutoff() : ModelParam() {}

  void set_model(LennardJonesAlpha * model) {
    model_ = model; }

  double compute(const int type1, const int type2,
      const ModelParams& model_params) override {
    const double cutoff = model_params.mixed_cutoff()[type1][type2];
    return model_->du_dr(cutoff, type1, type2, model_params);
  }

  EnergyDerivAtCutoff(std::istream& istr) : ModelParam(istr) {}
 private:
  /// temporary pointer for use with compute without requiring new arguments.
  LennardJonesAlpha * model_;
};

}  // namespace feasst

#endif  // FEASST_MODELS_LENNARD_JONES_ALPHA_H_
