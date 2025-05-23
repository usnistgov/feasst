
#ifndef FEASST_FLAT_HISTOGRAM_BIAS_H_
#define FEASST_FLAT_HISTOGRAM_BIAS_H_

#include <memory>
#include <string>
#include <map>

namespace feasst {

class CollectionMatrix;
class Histogram;
class LnProbability;
class Macrostate;

typedef std::map<std::string, std::string> argtype;

/**
  Bias for flat histogram Monte Carlo.
  Assumes a one-dimensional macrostate.

  CriteriaWriter outputs the following as a row for each Macrostate:
  - "ln_prob" the LnProbability.
  - "delta_ln_prob" the change in LnProbability from the previous Macrostate
    (e.g., ln_prob(state) - ln_prob(state - 1), so the lowest state is NaN).
 */
class Bias {
 public:
  Bias() {}

  /// Return the natural log of the bias for a transition from a macrostate
  /// in the old bin to a new bin.
  double ln_bias(const int bin_new, const int bin_old) const;

  /// Update only.
  virtual void update(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool is_endpoint,
    const Macrostate& macro) = 0;

  /// Perform an infrequent update to the bias.
  virtual void infrequent_update(const Macrostate& macro) {}

  /// The natural log of the macrostate probability.
  virtual const LnProbability& ln_prob() const = 0;

  virtual void resize(const Histogram& histogram) = 0;

  /// Set the macrostate probability distribution.
  virtual void set_ln_prob(const LnProbability& ln_prob) = 0;

  virtual std::string write() const { return std::string(""); }

  virtual std::string write_per_bin(const int bin) const;

  virtual std::string write_per_bin_header(const std::string& append) const;

  /// Return the simulation phase index used to differentiate production
  /// and initialization, etc.
  virtual int phase() const { return phase_; }

  /// Increment the simulation phase.
  virtual void increment_phase() { ++phase_; }

  /// Return the number of cycles required for completion.
  /// In TransitionMatrix and WLTM, this is the minimum number of sweeps.
  /// In WangLandau, this is the minimum number of flatness checks.
  /// Afterward, check again for completeness.
  virtual int cycles_to_complete() const = 0;

  /// Set the number of cycles required for completion.
  virtual void set_cycles_to_complete(const int cycle) = 0;

  virtual int num_cycles(const int state, const Macrostate& macro)
    const = 0;
  bool is_complete() const { return is_complete_; }
  void set_complete_() { is_complete_ = true; }

  // HWH hackish interface. See CollectionMatrixSplice::adjust_bounds.
  virtual void set_cm(const int macro, const Bias& bias);
  virtual const CollectionMatrix& cm() const;
  virtual const int visits(const int macro, const int index) const;
  virtual bool is_adjust_allowed(const Macrostate& macro) const {
    return false; }

  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Bias> create(std::istream& istr) const;
  virtual std::shared_ptr<Bias> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<Bias> >& deserialize_map();
  std::shared_ptr<Bias> deserialize(std::istream& istr);
  std::shared_ptr<Bias> factory(const std::string name, argtype * args);
  virtual ~Bias() {}

 protected:
  std::string class_name_ = "Bias";

  void set_incomplete_() { is_complete_ = false; }

  int bin_(
    const int macrostate_old,
    const int macrostate_new,
    const bool is_accepted);

  void serialize_bias_(std::ostream& ostr) const;
  explicit Bias(std::istream& istr);

 private:
  bool is_complete_ = false;
  int phase_ = 0;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_H_
