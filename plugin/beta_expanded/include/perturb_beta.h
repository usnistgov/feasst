
#ifndef FEASST_BETA_EXPANDED_PERTURB_BETA_H_
#define FEASST_BETA_EXPANDED_PERTURB_BETA_H_

#include "monte_carlo/include/perturb.h"

namespace feasst {

/**
  Change the inverse temperature, \f$\beta\f$, of the system by a fixed amount,
  randomly up or down.
 */
class PerturbBeta : public Perturb {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - fixed_beta_change: input the fixed amount that beta changes.
      The choice to increase or decrease is randomly selected.
   */
  explicit PerturbBeta(argtype args = argtype());
  explicit PerturbBeta(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held = false,
    Acceptance * acceptance = NULL) override;

  /// Change inverse temperature, \f$\beta\f$.
  void change_beta(const double delta_beta, System * system);

  void revert(System * system) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbBeta(std::istream& istr);
  virtual ~PerturbBeta() {}

  //@}
 private:
  double fixed_beta_change_;

  // temporary
  double previous_beta_;
};

inline std::shared_ptr<PerturbBeta> MakePerturbBeta(argtype args = argtype()) {
  return std::make_shared<PerturbBeta>(args);
}

}  // namespace feasst

#endif  // FEASST_BETA_EXPANDED_PERTURB_BETA_H_
