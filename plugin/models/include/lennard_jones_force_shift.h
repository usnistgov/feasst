
#ifndef FEASST_MODELS_LENNARD_JONES_FORCE_SHIFT_H_
#define FEASST_MODELS_LENNARD_JONES_FORCE_SHIFT_H_

#include "models/include/lennard_jones_cut_shift.h"

namespace feasst {

/**
  The Lennard-Jones potential, \f$U_{LJ}\f$ is described in LennardJonesAlpha.
  This class implements the force shifted (FS) version which ensures both
  \f$U(r_c)=0\f$ and \f$\left.\frac{\partial U}{\partial r}\right|_{r=r_c}=0\f$.

  \f$ U_{LJ}^{FS}(r) = \left\{
    \begin{array}{lr}
      U_{LJ}(r) - U_{LJ}(r_c) - \left.\frac{\partial U}{\partial r}\right|_{r=r_c} (r - r_c) & : r < r_c \\
      0 & r \ge r_c
    \end{array}
  \right. \f$
 */
class LennardJonesForceShift : public LennardJonesAlpha {
 public:
  explicit LennardJonesForceShift(argtype args = argtype());
  explicit LennardJonesForceShift(argtype * args);

  // HWH - optimize better: sqrt is used twice for distance
  /// Precompute the shift factor for optimization, given existing model parameters.
  void precompute(const ModelParams& existing) override;

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<LennardJonesForceShift>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<LennardJonesForceShift>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit LennardJonesForceShift(std::istream& istr);
  virtual ~LennardJonesForceShift() {}

 private:
  EnergyAtCutOff shift_;
  EnergyDerivAtCutOff force_shift_;
  bool precomputed_ = false;
};

inline std::shared_ptr<LennardJonesForceShift> MakeLennardJonesForceShift(
    argtype args = argtype()) {
  return std::make_shared<LennardJonesForceShift>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_LENNARD_JONES_FORCE_SHIFT_H_
