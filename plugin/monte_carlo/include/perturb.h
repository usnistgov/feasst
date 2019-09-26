
#ifndef FEASST_MONTE_CARLO_PERTURB_H_
#define FEASST_MONTE_CARLO_PERTURB_H_

#include <string>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/tunable.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

/**
  Perturbations to the system include the following types:

    1. move (single particle, regrowth, cluster moves).
    2. add/remove (grand canonical).
    3. change type (swap, growth expanded).
    4. change parameter (temperature, criteria swap bound).

  Perturbations may be followed up by one of the two following options:

    1. Revert the perturbation. For example, if a particle is moved, then
       restore the positions of that particle to the original.

    2. Finalize the perturbation. For example, if a selected particle is
       accepted for deletion, the particle is deleted in the finalize step.
 */
class Perturb {
 public:
  explicit Perturb(
    /**
      tunable_param : initial value of the tunable parameter (default: 0.1).
     */
    const argtype& args = argtype());

  /// Return the tunable parameter.
  const Tunable& tunable() const { return tunable_; }

  /// Set the minimum and maximum values of the tunable parameter.
  void set_tunable_min_and_max(const double min, const double max) {
    tunable_.set_min_and_max(min, max); }

  /// Set the value of the tunable parameter.
  void set_tunable(const double value) { tunable_.set_value(value); }

  /// Tune the parameter based on difference between target and actual.
  void tune(const double actual) { tunable_.tune(actual); }

  virtual void precompute(TrialSelect * select, System * system) {}

  /// Before perturbation, initialize some optimiation parameters.
  void before_select();

  /// Perturb the selection of the system.
  virtual void perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    /// If position is held, all but perform the actual perturbation, as
    /// typically required for calculations of old configurations and Rosenbluth
    const bool is_position_held = false
    ) { ERROR("not implemented"); }

  /// Return whether it is possible to revert.
  bool revert_possible() const { return revert_possible_; }

  /// Set whether it is possible to revert.
  void set_revert_possible(const bool revert_possible,
    /// If possible, store the selection.
    TrialSelect * revert_select);

  /// Return the stored revert selection.
  const TrialSelect * revert_select() const { return revert_select_; }

  /// Revert the system for an unaccepted perturbation.
  virtual void revert(System * system);

  /// Return whether it is possible to finalize.
  bool finalize_possible() const { return finalize_possible_; }

  /// Set whether it is possible to finalize.
  void set_finalize_possible(const bool finalize_possible,
    /// If possible, store the selection.
    TrialSelect * finalize_select);

  /// Return the stored finalize selection.
  const TrialSelect * finalize_select() const { return finalize_select_; }

  /// Finalize the system for an accepted perturbation.
  virtual void finalize(System * system);

  /// Return the probability.
  double probability() const { return probability_; }

  /// Set the probability.
  void set_probability(const double prob = 1) { probability_ = prob; }

  /// Print status header.
  virtual std::string status_header() const;

  /// Print status
  virtual std::string status() const;

  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Perturb> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Perturb> >& deserialize_map();
  std::shared_ptr<Perturb> deserialize(std::istream& istr);
  virtual ~Perturb() {}

 protected:
  std::string class_name_ = "Perturb";
  void disable_tunable_() { tunable_.disable(); }
  void serialize_perturb_(std::ostream& ostr) const;
  Perturb(std::istream& istr);

 private:
  Tunable tunable_;

  // optimzation or temporary object
  bool revert_possible_, finalize_possible_;
  TrialSelect * revert_select_;
  TrialSelect * finalize_select_;
  double probability_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_H_
