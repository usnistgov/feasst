
#ifndef FEASST_MONTE_CARLO_PERTURB_VOLUME_H_
#define FEASST_MONTE_CARLO_PERTURB_VOLUME_H_

#include <memory>
#include "monte_carlo/include/perturb.h"

namespace feasst {

class Select;

/**
  Change the volume of the system uniformly randomly in \f$\ln V\f$.
 */
class PerturbVolume : public Perturb {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - uniform_volume: if true, change volume uniformly in V instead of
      \f$\ln V\f$ (default: false).
    - constrain_volume_change: if true, use the previous volume change to
      do the opposite for use as the second stage in Gibbs ensemble
      (default: false).
    - Tunable arguments.
   */
  explicit PerturbVolume(argtype args = argtype());
  explicit PerturbVolume(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return true if the volume is changed uniformly.
  const bool uniform_volume() const { return uniform_volume_; }

  void perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held = false,
    Acceptance * acceptance = NULL) override;

  /// Change volume
  void change_volume(const double delta_volume,
      System * system,
      const Select& select);

  void precompute(TrialSelect * select, System * system) override;
  void revert(System * system) override;
  void finalize(System * system) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbVolume(std::istream& istr);
  virtual ~PerturbVolume() {}

  //@}
 private:
  bool uniform_volume_;
  bool constrain_volume_change_;
  argtype args_;

  // temporary
  double volume_change_;
};

inline std::shared_ptr<PerturbVolume> MakePerturbVolume(
    argtype args = argtype()) {
  return std::make_shared<PerturbVolume>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_VOLUME_H_
