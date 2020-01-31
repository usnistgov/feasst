
#ifndef FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_
#define FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "flat_histogram/include/bias.h"

namespace feasst {

/**
  Wang Landau flat histogram bias.
  https://doi.org/10.1103/PhysRevLett.86.2050
  https://doi.org/10.1063/1.1615966
 */
// HWH: Implement "gentle" WL where the bias is updated infrequently.
class WangLandau : public Bias {
 public:
  WangLandau(
    /**
      min_flatness : Number of flatness checks required for completion.

      add_to_ln_probability : The initial amount to add to the natural log of
        the macrostate probability upon visiting that state (default: 1.0).

      reduce_ln_probability : Reduce the amount to add to the natural log of the
        macrostate probability by multiplcation of this factor upon reaching a
        sufficiently flat histogram (default: 0.5).

      flatness_threshold : The visited states histogram is determined to be flat
        when the percentage between minimum visisted states and average reaches
        this threshold (default: 0.8).
     */
    const argtype &args = argtype());
  void update_or_revert(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool revert) override;
  const LnProbability& ln_prob() const override {
    return ln_prob_; }
  void resize(const Histogram& histogram) override;
  std::string write() const override;
  std::string write_per_bin(const int bin) const override;
  std::string write_per_bin_header() const override;
  void set_ln_prob(const LnProbability& ln_prob) override;
  std::shared_ptr<Bias> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  WangLandau(std::istream& istr);
  virtual ~WangLandau() {}

 private:
  LnProbability ln_prob_;
  double add_to_ln_probability_ = 0;
  double reduce_ln_probability_ = 0;
  double flatness_threshold_ = 0;

  /// Count of the number of times a state has been visited since the last time
  /// this histogram was reset after it was deemed to be sufficiently flat.
  std::vector<int> visited_states_;

  /// Number of times the visited states histogram was found to be flat.
  int num_flatness_ = 0;
  int min_flatness_ = 0;

  Arguments args_;

  /// Perform update when the visited states histogram is found to be flat.
  void flatness_update_();
};

inline std::shared_ptr<WangLandau> MakeWangLandau(
    const argtype& args = argtype()) {
  return std::make_shared<WangLandau>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_
