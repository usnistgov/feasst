
#ifndef FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_
#define FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "monte_carlo/include/criteria.h"

namespace feasst {

class Acceptance;
class Bias;
class LnProbability;
class Macrostate;
class Random;

typedef std::map<std::string, std::string> argtype;

/**
  Flat histogram acceptance criteria uses a bias to improve sampling and attempt
  to recover the free energy of the system as a function of the give macrostate.

  The Macrostate must be defined before the bias.

  CriteriaWriter outputs the following:
  - Criteria base class output.
  - Bias derived class output.
  - soft_min and soft_max, the current bounds of the macrostate.
  - rows of each Macrostate, given by "state," followed by Bias derived class-
    specific information for each Macrostate.
 */
class FlatHistogram : public Criteria {
 public:
  FlatHistogram();  // do not use this constructor.

  //@{
  /** @name Arguments
    - Macrostate: MacrostateNumParticles, MacrostateEnergy, etc
    - Bias: WangLandau, TransitionMatrix, WLTM, etc.
    - Criteria arguments.
   */
  explicit FlatHistogram(argtype args);
  explicit FlatHistogram(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Construct with a macrostate and a bias
  FlatHistogram(std::shared_ptr<Macrostate> macrostate,
      std::shared_ptr<Bias> bias);

  /// Same as above, but with an added Constraint.
  FlatHistogram(std::shared_ptr<Macrostate> macrostate,
      std::shared_ptr<Bias> bias,
      std::shared_ptr<Constraint> constraint);

  /// Return the macrostate.
  const Macrostate& macrostate() const override;

  Macrostate * get_macrostate() { return macrostate_.get(); }

  /// Return the bias.
  const Bias& bias() const override;
  void set_bias(std::shared_ptr<Bias> bias);

  int cycles_to_complete() const override;
  void set_cycles_to_complete(const int num) override;
  int num_cycles(const int state = -1) const override;
  bool is_complete() const override;
  void set_complete() override;

  void before_attempt(const System& system) override;

  bool is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) override;

  std::string write() const override;
  int phase() const override;
  void increment_phase() override;

  /// Return the state. Return -1 if state is not determined.
  int state() const override;
  int num_states() const override;
  int state_old() const override;
  int state_new() const override;
  void update_state(const System& system, const Acceptance& accept) override;

  /// Set the macrostate probability distribution.
  void set_ln_prob(const LnProbability& ln_prob);

  /// Return the macrostate probability distribution.
  const LnProbability& ln_prob() const;

  // HWH hackish implementation for prefetch
  // Revert changes from previous trial.
  void revert_(const bool accepted, const bool endpoint,
               const double ln_prob, const std::vector<int>& updated) override;
  // HWH rename: delete
  void finalize(const Acceptance& acceptance) override;
  void revert(const Acceptance& acceptance) override { finalize(acceptance); }
  void imitate_trial_rejection_(const double ln_prob,
      const int state_old,
      const int state_new,
      const bool endpoint) override;
  void update() override;
  bool is_fh_equal(const FlatHistogram& flat_histogram,
    const double tolerance) const;

  // HWH hackish adjust_bounds interface. See CollectionMatrixSplice.
  int set_soft_max(const int index, const System& sys) override;
  int set_soft_min(const int index, const System& sys) override;
  void set_cm(const bool inc_max, const int macro,
              const Criteria& crit) override;
  void adjust_bounds(const bool left_most, const bool right_most,
    const bool left_complete, const bool right_complete,
    const bool all_min_size, const int min_size, const System& system,
    const System * upper_sys, Criteria * criteria, bool * adjusted_up,
    std::vector<int> * states) override;
  const FlatHistogram& flat_histogram() const override { return *this; }
  int soft_min() const override;
  int soft_max() const override;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<FlatHistogram>(istr); }
  std::shared_ptr<Criteria> create(argtype * args) const override {
    return std::make_shared<FlatHistogram>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit FlatHistogram(std::istream& istr);
  std::unique_ptr<FlatHistogram> flat_histogram(const Criteria& criteria);
  ~FlatHistogram();

  //@}
 private:
  std::shared_ptr<Bias> bias_;
  std::shared_ptr<Macrostate> macrostate_;
  int macrostate_old_ = -1;
  int macrostate_new_ = -1;
  int macrostate_current_ = -1;
  bool is_macrostate_set_ = false;

  // temporary
  std::unique_ptr<Acceptance> empty_;

  void init_(std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias);
  void check_left_and_right_most_(const bool left_most, const bool right_most,
    const bool all_min_size,
    const int min_size, const System& system, const System * upper_sys,
    Criteria * criteria);
};

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram() {
  return std::make_shared<FlatHistogram>();
}

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram(
    std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias) {
  return std::make_shared<FlatHistogram>(macrostate, bias);
}

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram(argtype args) {
  return std::make_shared<FlatHistogram>(args);
}

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram(
    std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias,
    std::shared_ptr<Constraint> constraint) {
  return std::make_shared<FlatHistogram>(macrostate, bias, constraint);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_
