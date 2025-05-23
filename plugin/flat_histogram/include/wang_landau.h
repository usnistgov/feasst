
#ifndef FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_
#define FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_

#include <memory>
#include <string>
#include <vector>
#include "flat_histogram/include/bias.h"

namespace feasst {

class LnProbability;

/**
\rst
Wang-Landau flat histogram bias.\ :footcite:p:`wang_efficient_2001`
\endrst

  CriteriaWriter outputs the following:
  - num_flatness: The number of times the visited states histogram was found
    to be flat.
  - rows for each Macrostate with the following:
    - "visited" the number of times a state has been visited since the last
      time this visited states histogram was reset after it was deemed to be
      sufficiently flat.

\rst
References:

.. footbibliography::
\endrst
 */
class WangLandau : public Bias {
 public:
  //@{
  /** @name Arguments
    - min_flatness : Number of flatness checks required for completion.
    - add_to_ln_probability : The initial amount to add to the natural log of
      the macrostate probability upon visiting that state (default: 1.0).
    - reduce_ln_probability : Reduce the amount to add to the natural log of the
      macrostate probability by multiplcation of this factor upon reaching a
      sufficiently flat histogram (default: 0.5).
    - flatness_threshold : The visited states histogram is determined to be flat
      when the percentage between minimum visisted states and average reaches
      this threshold (default: 0.8).
    - min_visit_per_macro: The minimum number of visits for each macrostate
      required during flatness check (default: 10^3).
   */
  explicit WangLandau(argtype args = argtype());
  explicit WangLandau(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  int min_flatness() const { return min_flatness_; }
  void update(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool is_endpoint,
    const Macrostate& macro) override;
  int cycles_to_complete() const override { return min_flatness_; }
  void set_cycles_to_complete(const int flatness) override;
  int num_cycles(const int state, const Macrostate& macro) const override {
    return num_flatness_;}
  const LnProbability& ln_prob() const override;
  void resize(const Histogram& histogram) override;
  void infrequent_update(const Macrostate& macro) override;
  std::string write() const override;
  std::string write_per_bin(const int bin) const override;
  std::string write_per_bin_header(const std::string& append) const override;
  void set_ln_prob(const LnProbability& ln_prob) override;
  const int num_flatness() const { return num_flatness_; }
  std::shared_ptr<Bias> create(std::istream& istr) const override {
    return std::make_shared<WangLandau>(istr); }
  std::shared_ptr<Bias> create(argtype * args) const override {
    return std::make_shared<WangLandau>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit WangLandau(std::istream& istr);
  virtual ~WangLandau();

  //@}
 private:
  std::unique_ptr<LnProbability> ln_prob_;
  double add_to_ln_probability_ = 0;
  double reduce_ln_probability_ = 0;
  double flatness_threshold_ = 0;
  int min_visit_per_macro_;

  std::vector<int> visited_states_;

  /// Number of times the visited states histogram was found to be flat.
  int num_flatness_ = 0;
  int min_flatness_ = 0;

  void flatness_check_();
  /// Perform update when the visited states histogram is found to be flat.
  void flatness_update_();
};

inline std::shared_ptr<WangLandau> MakeWangLandau(argtype args = argtype()) {
  return std::make_shared<WangLandau>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_
