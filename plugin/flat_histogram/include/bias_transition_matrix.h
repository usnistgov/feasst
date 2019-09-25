
#ifndef FEASST_FLAT_HISTOGRAM_BIAS_TRANSITION_MATRIX_H_
#define FEASST_FLAT_HISTOGRAM_BIAS_TRANSITION_MATRIX_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "flat_histogram/include/bias.h"
#include "flat_histogram/include/collection_matrix.h"

namespace feasst {

/**
  Transition matrix flat histogram bias.

  Direct calculation of liquid–vapor phase equilibria from transition matrix
  Monte Carlo simulation
  https://doi.org/10.1063/1.1572463

  Elucidating the effects of adsorbent flexibility on fluid adsorption using
  simple models and flat-histogram sampling methods.
  https://doi.org/10.1063/1.4884124
 */
class BiasTransitionMatrix : public Bias {
 public:
  BiasTransitionMatrix(
    /**
      min_visits : A sweep is performed when all macrostates are visited by
        another macrostate this number of times (default: 100).

      min_sweeps : Number of sweeps required for completion.

      num_steps_to_update : Number of steps to update macrostate
        distribution (default: 1e6).
     */
    const argtype &args = argtype());
  void update(const int macrostate_old,
              const int macrostate_new,
              const double ln_metropolis_prob,
              const bool is_accepted) override;
  const LnProbabilityDistribution& ln_macro_prob() const override {
    return ln_macro_prob_; }
  void resize(const Histogram& histogram) override;
  void revert(const int macrostate_new, const int macrostate_old) override;
  std::string write() const override;
  std::string write_per_bin(const int bin) const override;
  std::string write_per_bin_header() const override;

  std::shared_ptr<Bias> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit BiasTransitionMatrix(std::istream& istr);
  virtual ~BiasTransitionMatrix() {}

 protected:
  void infrequent_update_();

 private:
  std::string class_name_ = "BiasTransitionMatrix";
  LnProbabilityDistribution ln_macro_prob_;
  TripleBandedCollectionMatrix collection_;
  std::vector<int> visits_;
  int min_visits_ = 0;
  int num_sweeps_ = 0;
  int min_sweeps_ = 0;
  int num_steps_to_update_ = 0;
  int num_steps_since_update_ = 0;

  std::vector<BiasTransitionMatrix> blocks_;
  bool is_block_ = false;
  int iter_block_= -1;

  void update_blocks_(
      const int macrostate_old,
      const int macrostate_new,
      const double ln_metropolis_prob,
      const bool is_accepted);
};

inline std::shared_ptr<BiasTransitionMatrix> MakeBiasTransitionMatrix(
    const argtype &args = argtype()) {
  return std::make_shared<BiasTransitionMatrix>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_TRANSITION_MATRIX_H_
