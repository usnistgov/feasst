
#ifndef FEASST_FLAT_HISTOGRAM_TRANSITION_MATRIX_GUESS_H_
#define FEASST_FLAT_HISTOGRAM_TRANSITION_MATRIX_GUESS_H_

#include <memory>
#include <string>
#include <vector>
#include "flat_histogram/include/bias.h"

namespace feasst {

class LnProbability;
class TransitionMatrix;

/**
  Begin with a guess of the LnProbability Macrostate distribution, and end with
  TransitionMatrix.
*/
class TransitionMatrixGuess : public Bias {
 public:
  //@{
  /** @name Arguments
    - TransitionMatrix arguments.
    - min_collect_sweeps: Begin using the TransitionMatrix bias when this many
      sweeps has been reached (default: 10).
    - ln_prob_file: name of file containing one number for each line
      corresponding to the LnProbability of each Macrostate.
   */
  explicit TransitionMatrixGuess(argtype args = argtype());
  explicit TransitionMatrixGuess(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void update(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool is_endpoint,
    const Macrostate& macro) override;

  /// Updates min_sweeps, but neither flatness.
  int cycles_to_complete() const override;
  void set_cycles_to_complete(const int sweeps) override;
  int num_cycles(const int state, const Macrostate& macro) const override;
  const TransitionMatrix& transition_matrix() const;
  const LnProbability& ln_prob() const override;
  void resize(const Histogram& histogram) override;
  void infrequent_update(const Macrostate& macro) override;
  std::string write() const override;
  std::string write_per_bin(const int bin) const override;
  std::string write_per_bin_header(const std::string& append) const override;
  void set_ln_prob(const LnProbability& ln_prob) override;

  // HWH hackish interface. See CollectionMatrixSplice::adjust_bounds.
  void set_cm(const int macro, const Bias& bias) override;
  const CollectionMatrix& cm() const override;
  const int visits(const int macro, const int index) const override;
  bool is_adjust_allowed(const Macrostate& macro) const override;

  std::shared_ptr<Bias> create(std::istream& istr) const override;
  std::shared_ptr<Bias> create(argtype * args) const override {
    return std::make_shared<TransitionMatrixGuess>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TransitionMatrixGuess(std::istream& istr);
  virtual ~TransitionMatrixGuess();

  //@}
 private:
  int min_collect_sweeps_;
  int production_ = 0;
  std::unique_ptr<LnProbability> ln_prob_guess_;
  std::unique_ptr<TransitionMatrix> transition_matrix_;
  bool is_tm_bias_(const Macrostate& macro) const;

  // temporary and not serialized
  bool is_tm_bias_at_update_ = false;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_TRANSITION_MATRIX_GUESS_H_
