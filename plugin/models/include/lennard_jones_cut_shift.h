
#ifndef FEASST_MODELS_LENNARD_JONES_CUT_SHIFT_H_
#define FEASST_MODELS_LENNARD_JONES_CUT_SHIFT_H_

#include "models/include/lennard_jones_alpha.h"

namespace feasst {

/**
  The Lennard-Jones potential, \f$U_{LJ}\f$ is described in LennardJonesAlpha.
  This class implements the cut and shifted (CS) version which ensures \f$U(r_c)=0\f$.

  \f$ U_{LJ}^{CS}(r) = \left\{
    \begin{array}{lr}
      U_{LJ}(r) - U_{LJ}(r_c) & : r < r_c \\
      0 & : r \ge r_c
    \end{array}
  \right. \f$

 For a Weeks-Chandler-Anderson (WCA) potential, use this class.
 The cutoffs are computed by LennardJonesCutShift::set_wca.
 Thus, one workflow is as follows:

   auto wca = MakeLennardJonesCutShift();
   ModelParams wca_params = configuration.model_params(); // copy model params
   wca->set_wca(type1, type2, wca_params);  // set cutoff
   wca->precompute(wca_params);   // compute shifts, etc
   auto potential = MakePotential(wca);
   potential->set(wca_params); // use wca_params.

 */
class LennardJonesCutShift : public LennardJonesAlpha {
 public:
  explicit LennardJonesCutShift(argtype args = argtype());
  explicit LennardJonesCutShift(argtype * args);

  // HWH some issues with this implementation include
  // - what if model params change, or is defined different by a special potential
  // - how to simplify user interface
  /// Precompute the shift factor for optimization, given existing model parameters.
  void precompute(const Configuration& config) override;

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<LennardJonesCutShift>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<LennardJonesCutShift>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit LennardJonesCutShift(std::istream& istr);
  virtual ~LennardJonesCutShift() {}

private:
  EnergyAtCutOff shift_;
};

inline std::shared_ptr<LennardJonesCutShift> MakeLennardJonesCutShift(
    argtype args = argtype()) {
  return std::make_shared<LennardJonesCutShift>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_LENNARD_JONES_CUT_SHIFT_H_
