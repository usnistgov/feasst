
#ifndef FEASST_FLAT_HISTOGRAM_TRANSITION_MATRIX_H_
#define FEASST_FLAT_HISTOGRAM_TRANSITION_MATRIX_H_

#include <memory>
#include <string>
#include <vector>
#include "flat_histogram/include/bias.h"

namespace feasst {

class CollectionMatrix;
class Histogram;
class LnProbability;

/**
\rst
Transition matrix\ :footcite:p:`errington_direct_2003` flat histogram bias
may be used to compute phase behavior in the grand canonical ensemble.

Elucidating the effects of adsorbent flexibility on fluid adsorption using
simple models and flat-histogram sampling methods.
https://doi.org/10.1063/1.4884124

A count of all accepted transitions between macrostates is used for a newly
developed sweep metric.

Note that transitions to macrostates that have not yet been sampled are
unbiased.
This can lead to long initialization times, depending on your choice of the
chemical potential.
For example, if a simulation begins with zero particles and attempts a
particle insertion, this insertion is unlikely to be accepted if the chemical
potential is very low (favoring less particles).
This can lead to very low initialization times.
This issue can be overcome with the use of WLTM, which instead uses
Wang-Landau to initialize the transition matrix.
This makes WLTM efficiency less sensitive than TransitionMatrix to the choice
of chemical potential.

References:

.. footbibliography::
\endrst
 */
class TransitionMatrix : public Bias {
 public:
  //@{
  /** @name Arguments
    - CollectionMatrix arguments.
    - min_visits: A sweep is performed when all macrostates are visited by
      another macrostate this number of times (default: 100).
    - average_visits: A sweep is performed when macrostates are visited by
      another macrostate more than this average number of times (default: 0).
    - min_sweeps: Number of sweeps required for completion.
    - reset_sweeps: The 'phase' counter increments from 0 to 1 when
      reset_sweeps are completed (default: -1 [counter will never increment])
    - new_sweep: if set to 1, use new sweep definition of "the minimum number
      of accepted transitions for each possible" (default: 0).
    - widom: if true, also store a CollectionMatrix without the minimum(1,P)
      term typical of Metropolis acceptance (default: false).
      This allows direct calculation of the chemical potential
   */
  explicit TransitionMatrix(argtype args = argtype());
  explicit TransitionMatrix(argtype * args);

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

  /// Return the minimum sweeps required for completion.
  int min_sweeps() const { return min_sweeps_; }

  /// Return the reset_sweeps parameter.
  int reset_sweeps() const { return reset_sweeps_; }

  int cycles_to_complete() const override { return min_sweeps_;}
  void set_cycles_to_complete(const int sweeps) override;
  int num_cycles(const int state, const Macrostate& macro) const override;
  const LnProbability& ln_prob() const override;
  void resize(const int size);
  void resize(const Histogram& histogram) override;
  std::string write() const override;
  std::string write_per_bin(const int bin) const override;
  std::string write_per_bin_header() const override;
  void set_ln_prob(const LnProbability& ln_prob) override;
  void infrequent_update(const Macrostate& macro) override;
  bool is_equal(const TransitionMatrix& transition_matrix,
    const double tolerance) const;

  // HWH hackish interface. See CollectionMatrixSplice::adjust_bounds.
  void set_cm(const CollectionMatrix& cm);
  void set_cm(const int macro, const Bias& bias) override;
  const CollectionMatrix& cm() const override;
  const int visits(const int macro, const int index) const override;
  bool is_adjust_allowed(const Macrostate& macro) const override;

  std::shared_ptr<Bias> create(std::istream& istr) const override;
  std::shared_ptr<Bias> create(argtype * args) const override {
    return std::make_shared<TransitionMatrix>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TransitionMatrix(std::istream& istr);
  // std::unique_ptr<TransitionMatrix> tm(const Bias& bias);
  virtual ~TransitionMatrix();

  //@}
 private:
  std::unique_ptr<CollectionMatrix> collection_;
  std::unique_ptr<CollectionMatrix> collection_widom_;
  std::unique_ptr<LnProbability> ln_prob_;
  std::vector<std::vector<int> > visits_;
  int min_visits_ = 0;
  int num_sweeps_ = 0;
  int min_sweeps_ = 0;
  int reset_sweeps_ = -1;
  int average_visits_ = 0;
  int new_sweep_ = 0;
  bool widom_;

  // find minimum visits in soft range
  int min_vis_calc_(const Macrostate& macro) const;
};

inline std::shared_ptr<TransitionMatrix> MakeTransitionMatrix(
    argtype args = argtype()) {
  return std::make_shared<TransitionMatrix>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_TRANSITION_MATRIX_H_
