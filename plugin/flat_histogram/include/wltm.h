
#ifndef FEASST_FLAT_HISTOGRAM_WLTM_H_
#define FEASST_FLAT_HISTOGRAM_WLTM_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/transition_matrix.h"

namespace feasst {

/**
  Begin with WangLandau and end with TransitionMatrix.
 */
class WLTM : public Bias {
 public:
  /**
    args:
    - collect_flatness: Begin populating the collection matrix when Wang-Landau
      has completed this many flatness checks.
      Note that populating the collection matrix does not necessarily mean that
      the collection matrix is used to compute the bias.
    - min_flatness: After this many flatness checks have been completed
      with Wang-Landau, use the bias from the collection matrix instead.
      Also, increment the phase when this occurs.
    - min_sweeps: Number of sweeps required for completion.
   */
  WLTM(const argtype &args = argtype());
  void update_or_revert(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool revert) override;
  void set_num_iterations(const int sweeps) override {
    transition_matrix_->set_num_iterations(sweeps); }
  const LnProbability& ln_prob() const override;
  void resize(const Histogram& histogram) override;
  void infrequent_update() override;
  std::string write() const override;
  std::string write_per_bin(const int bin) const override;
  std::string write_per_bin_header() const override;
  void set_ln_prob(const LnProbability& ln_prob) override;
  std::shared_ptr<Bias> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit WLTM(std::istream& istr);
  virtual ~WLTM() {}

 private:
  int collect_flatness_;
  int min_flatness_;
  int production_ = 0;
  std::shared_ptr<WangLandau> wang_landau_;
  std::shared_ptr<TransitionMatrix> transition_matrix_;
};

inline std::shared_ptr<WLTM> MakeWLTM(
    const argtype& args = argtype()) {
  return std::make_shared<WLTM>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WLTM_H_
