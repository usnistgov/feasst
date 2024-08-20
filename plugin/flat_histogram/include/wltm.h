
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
  Begin with WangLandau and end with TransitionMatrix.
  See https://doi.org/10.1063/1.1615966 and https://doi.org/10.1063/1.4884124
 */
class WLTM : public Bias {
 public:
  //@{
  /** @name Arguments
    - WangLandau arguments.
    - TransitionMatrix arguments.
    - collect_flatness: Begin populating the collection matrix when Wang-Landau
      has completed this many flatness checks.
      Note that populating the collection matrix does not necessarily mean that
      the collection matrix is used to compute the bias.
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
  int num_iterations_to_complete() const override;
  void set_num_iterations_to_complete(const int sweeps) override;
  int num_iterations(const int state, const Macrostate& macro) const override;
  const TransitionMatrix& transition_matrix() const;
  const LnProbability& ln_prob() const override;
  void resize(const Histogram& histogram) override;
  void infrequent_update(const Macrostate& macro) override;
  std::string write() const override;
  std::string write_per_bin(const int bin) const override;
  std::string write_per_bin_header() const override;
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
};

inline std::shared_ptr<WLTM> MakeWLTM(argtype args = argtype()) {
  return std::make_shared<WLTM>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WLTM_H_
