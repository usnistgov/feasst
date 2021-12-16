
#ifndef FEASST_MODELS_LENNARD_JONES_ALPHA_H_
#define FEASST_MODELS_LENNARD_JONES_ALPHA_H_

#include "utils/include/arguments.h"
#include "configuration/include/model_params.h"
#include "system/include/lennard_jones.h"

namespace feasst {

/**
  This Lennard-Jones potential, \f$U_{LJ}\f$ is a generalization of the one
  described in LennardJones.

  \f$ U_{LJ} = 4\epsilon \left[ \left(\frac{\sigma + \Delta_\sigma}{r + \Delta_\sigma}\right)^{2\alpha}
                              - \left(\frac{\sigma + \Delta_\sigma}{r + \Delta_\sigma}\right)^\alpha \right] \f$,

  In this class, the \f$\alpha\f$ parameter may be varied (default: 6).
  In addition, the delta_sigma parameter may be set as a site-specific
  ModelParam in the fstprt file (default: 0).

  Note that \f$U_{LJ}\f$ crosses zero when \f$r=\sigma\f$.

  Consider an alternative interpretation for \f$\Delta_\sigma = \sigma_r - \sigma\f$

  \f$ U_{LJ} = 4\epsilon \left[ \left(\frac{\sigma_r}{r - \sigma + \sigma_r}\right)^{2\alpha}
                              - \left(\frac{\sigma_r}{r - \sigma + \sigma_r}\right)^\alpha \right] \f$

  Thus, \f$\sigma_r = \sigma + \Delta_\sigma\f$ determines the well width and shape of the potential,
  while \f$\sigma\f$ determines the excluded volume by shifting the potential right or left.

  An additional option is the use of a "lambda" potential which decouples the
  attractive and repulsive portions of the LJ potential to allow for a modified
  well-depth or even a shoulder, as shown in
  Fig. 1 of https://doi.org/10.1039/C7SM01005B
  and also described in https://dx.doi.org/10.1021/ja802124e.

  \f$ U^\lambda = \left\{
    \begin{array}{lr}
      U_{LJ} + \epsilon(1 - \lambda_{ij}) & r \le r_m \\
      \lambda_{ij}U_{LJ} & r > r_m
    \end{array}
  \right. \f$

  where \f$r_m = (\sigma + \Delta_\sigma)2^{1/\alpha}-\Delta_\sigma\f$ is the minimum.
  The lambda parameter is also a site-specific ModelParam in the fstprt file.
 */
class LennardJonesAlpha : public LennardJones {
 public:
  /**
    args:
    - alpha: set the value of \f$\alpha\f$ (default: 6).
    - lambda: enable "lambda" potentials (default: false)
   */
  explicit LennardJonesAlpha(argtype args = argtype());
  explicit LennardJonesAlpha(argtype * args);

  /// Return the value of alpha.
  const double& alpha() const { return alpha_; }

  /// Initialize WCA cutoff distances for types 1 and 2 in model parameters.
  void set_wca(const int site_type1, const int site_type2,
      ModelParams * params) const;

  void precompute(const ModelParams& existing) override;

  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) override;

  /// Return the derivative in the potential energy with respect to distance.
  double du_dr(
    const double distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) const;

  // Derived classes may call the LennardJonesAlpha energy directly.
  double energy_without_shift(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) {
    return LennardJonesAlpha::energy(squared_distance, type1, type2, model_params);
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<LennardJonesAlpha>(istr); }

  void serialize(std::ostream& ostr) const override;
  explicit LennardJonesAlpha(std::istream& istr);
  virtual ~LennardJonesAlpha() {}

 protected:
  void serialize_lennard_jones_alpha_(std::ostream& ostr) const;

 private:
  double alpha_;
  int delta_sigma_index_ = -1;
  bool lambda_;
  int lambda_index_ = -1;
  double two_raised_inv_alpha_;
};

inline std::shared_ptr<LennardJonesAlpha> MakeLennardJonesAlpha(
    argtype args = argtype()) {
  return std::make_shared<LennardJonesAlpha>(args);
}

class DeltaSigma : public ModelParam {
 public:
  DeltaSigma() : ModelParam() { set_name("delta_sigma"); }
  DeltaSigma(std::istream& istr) : ModelParam(istr) {}
};

inline std::shared_ptr<DeltaSigma> MakeDeltaSigma() {
  return std::make_shared<DeltaSigma>();
}

class Lambda : public ModelParam {
 public:
  Lambda() : ModelParam() { set_name("lambda"); }
  Lambda(std::istream& istr) : ModelParam(istr) {}
};

inline std::shared_ptr<Lambda> MakeLambda() {
  return std::make_shared<Lambda>();
}

class EnergyAtCutoff : public ModelParam {
 public:
  EnergyAtCutoff() : ModelParam() { set_name("energy_at_cutoff"); }

  void set_model(LennardJonesAlpha * model) {
    model_ = model; }

  double compute(const int type1, const int type2,
      const ModelParams& model_params) override;

  EnergyAtCutoff(std::istream& istr) : ModelParam(istr) {}

private:
  /// temporary pointer for use with compute without requiring new arguments.
  LennardJonesAlpha * model_;
};

class EnergyDerivAtCutoff : public ModelParam {
 public:
  EnergyDerivAtCutoff() : ModelParam() { set_name("energy_deriv_at_cutoff"); }

  void set_model(LennardJonesAlpha * model) {
    model_ = model; }

  double compute(const int type1, const int type2,
      const ModelParams& model_params) override;

  EnergyDerivAtCutoff(std::istream& istr) : ModelParam(istr) {}
 private:
  /// temporary pointer for use with compute without requiring new arguments.
  LennardJonesAlpha * model_;
};

}  // namespace feasst

#endif  // FEASST_MODELS_LENNARD_JONES_ALPHA_H_
