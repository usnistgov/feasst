
#ifndef FEASST_MAYER_MAYER_SAMPLING_H_
#define FEASST_MAYER_MAYER_SAMPLING_H_

#include "monte_carlo/include/criteria.h"
#include "math/include/accumulator.h"

namespace feasst {

class Random;

/**
  Mayer-sampling Monte Carlo acceptance criteria (see
  https://doi.org/10.1103/PhysRevLett.92.220601).
 */
class MayerSampling : public Criteria {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - num_trials_per_iteration: define an iteration as a number of trials
      (as measured by number of calls to is_accepted) default: 1e9.
    - intra_potential: index of intramolecular potential that will be used
      to select the move. Ignore if -1 (default: -1).
    - Criteria arguments.
   */
  explicit MayerSampling(argtype args = argtype());
  explicit MayerSampling(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Disable optimization when overlap is detected.
  void precompute(System * system) override;

  bool is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) override;

  /// Return the mayer ensemble of the full potential.
  const Accumulator& mayer() const { return mayer_; }

  /// Return the mayer ensemble of the reference potential.
  const Accumulator& mayer_ref() const { return mayer_ref_; }

  /// Return the ratio of the second virial coefficient of the full potential
  /// to the second virial coefficient of the reference potential.
  double second_virial_ratio() const;

  /**
    Return the block standard deviation of the second_virial_ratio.
    This is computed using the error propogation formula variance formula:

    \f$f = mayer/mayer_ref\f$

    \f$\sigma_f = \sqrt{(\sigma_{mayer}/mayer_ref)^2 + (f\sigma_{mayer_ref}/mayer_ref)^2}\f$
   */
  double second_virial_ratio_block_stdev() const;

  std::string write() const override;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<MayerSampling>(istr); }
  std::shared_ptr<Criteria> create(argtype * args) const override {
    return std::make_shared<MayerSampling>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit MayerSampling(std::istream& istr);
//  explicit MayerSampling(const Criteria& criteria);
  ~MayerSampling() {}

  //@}
 private:
  double f12old_ = -1.;
  double f12ref_ = -1.;
  int num_trials_per_iteration_;
  Accumulator mayer_;
  Accumulator mayer_ref_;
  int intra_pot_;
};

inline std::shared_ptr<MayerSampling> MakeMayerSampling(
    argtype args = argtype()) {
  return std::make_shared<MayerSampling>(args);
}

}  // namespace feasst

#endif  // FEASST_MAYER_MAYER_SAMPLING_H_
