
#ifndef FEASST_FLAT_HISTOGRAM_WLTM_H_
#define FEASST_FLAT_HISTOGRAM_WLTM_H_

#include <memory>
#include <string>
#include <vector>
#include "flat_histogram/include/bias.h"

namespace feasst {

class WangLandau;
class TransitionMatrix;

/**
\rst
Begin with WangLandau and end with TransitionMatrix.\ :footcite:p:`shell_improved_2003`
For the benefits of using WangLandau for initialization and
TransitionMatrix for production convergence, see the Appendix of :footcite:t:`shen_elucidating_2014`.
\endrst

  WLTM operates in the follow three stages.

  1. WangLandau only.
     With CriteriaWriter, WLTM writes both WangLandau and TransitionMatrix.
     At this stage, the LnProbability from WangLandau appears as ln_prob.
     Because TransitionMatrix has not been used yet, the bias appears as all
     zeros with the header of ln_prob_tm.

  2. When the WangLandau flatness reaches collect_flatness, the
     CollectrionMatrix begins to update with each trial move, but the bias is
     still based off of WangLandau.

  3. When the WangLandau flatness reaches min_flatness, and the number of
     CollectionMatrix sweeps is greater than min_collect_sweeps, then the bias
     switches to TransitionMatrix.
     Writh CriteriaWriter, the TransitionMatrix LnProbability appears as
     ln_prob, while the WangLandau bias appears as ln_prob_wl.

\rst
References:

.. footbibliography::
\endrst
*/
class WLTM : public Bias {
 public:
  //@{
  /** @name Arguments
    - WangLandau arguments.
    - TransitionMatrix arguments.
    - collect_flatness: Begin populating the CollectionMatrix when WangLandau
      has completed this many flatness checks.
      Note that populating the CollectionMatrix does not necessarily mean that
      the CollectionMatrix is used to compute the Bias.
    - min_collect_sweeps: In addition to WangLandau::min_flatness, do not use
      TransitionMatrix as bias until it has this minimum number of sweeps.
      If -1, do nothing (default: -1).
   */
  explicit WLTM(argtype args = argtype());
  explicit WLTM(argtype * args);

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
    return std::make_shared<WLTM>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit WLTM(std::istream& istr);
  virtual ~WLTM();

  //@}
 private:
  int collect_flatness_;
  int min_flatness_;
  int min_collect_sweeps_;
  int production_ = 0;
  std::unique_ptr<WangLandau> wang_landau_;
  std::unique_ptr<TransitionMatrix> transition_matrix_;

  bool is_wl_bias_(const Macrostate& macro) const;
  bool is_cm_update_() const;

  // temporary and not serialized
  bool is_wl_bias_at_update_ = false;
};

inline std::shared_ptr<WLTM> MakeWLTM(argtype args = argtype()) {
  return std::make_shared<WLTM>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WLTM_H_
