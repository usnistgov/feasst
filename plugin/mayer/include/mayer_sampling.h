
#ifndef FEASST_MAYER_MAYER_SAMPLING_H_
#define FEASST_MAYER_MAYER_SAMPLING_H_

#include "monte_carlo/include/criteria.h"
#include "math/include/accumulator.h"

namespace feasst {

class Random;

/**
  Mayer-sampling Monte Carlo acceptance criteria (see
  https://doi.org/10.1103/PhysRevLett.92.220601).

  Extrapolation in temperature is performed as descirbed in:
  https://doi.org/10.1063/1.5016165 .

  The coefficients for the Taylor series extrapolation are written to file
  and include the division by factorial.
 */
class MayerSampling : public Criteria {
 public:
  //@{
  /** @name Arguments
    - num_trials_per_iteration: define an iteration as a number of trials
      (as measured by number of calls to is_accepted) default: 1e9.
    - intra_potential: index of intramolecular potential that will be used
      to select the move. Ignore if -1 (default: -1).
    - num_beta_taylor: number of derivatives of second virial ratio with
      respect to beta. (default: 0).
    - training_file: if not empty, file name to write training data
      (default: empty).
    - training_per_write: write every this many sets of data (default: 1e4).
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

  /// Return the number of beta derivatives, starting with 1.
  int num_beta_taylor() const { return static_cast<int>(beta_taylor_.size()); }

  /// Return the beta derivatives in the mayer function f12
  const std::vector<Accumulator>& beta_taylor() const {
    return beta_taylor_; }

  /// Return the beta derivative in the second virial ratio by n factorial
  double beta_taylor(const int deriv) const;

  std::string write() const override;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<MayerSampling>(istr); }
  std::shared_ptr<Criteria> create(argtype * args) const override {
    return std::make_shared<MayerSampling>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit MayerSampling(std::istream& istr);
//  explicit MayerSampling(const Criteria& criteria);
  virtual ~MayerSampling() {}

  //@}
 private:
  double f12old_ = -1.;
  double f12ref_ = -1.;
  int num_trials_per_iteration_;
  Accumulator mayer_;
  Accumulator mayer_ref_;
  int intra_pot_;
  std::vector<Accumulator> beta_taylor_;

  // for training output
  std::string training_file_;
  std::vector<std::vector<double> > data_;
  int training_per_write_;
};

inline std::shared_ptr<MayerSampling> MakeMayerSampling(
    argtype args = argtype()) {
  return std::make_shared<MayerSampling>(args);
}

}  // namespace feasst

#endif  // FEASST_MAYER_MAYER_SAMPLING_H_
