
#ifndef FEASST_MONTE_CARLO_ANALYZE_H_
#define FEASST_MONTE_CARLO_ANALYZE_H_

#include <vector>
#include <memory>
#include <string>
#include <map>
#include "monte_carlo/include/stepper.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Perform a read-only action every so many steps.
 */
class Analyze : public Stepper {
 public:
  explicit Analyze(const argtype &args = argtype()) : Stepper(args) {}

  /// Initialize and precompute before trials.
  virtual void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) {}

  /// Check every trial if action is to be performed.
  virtual void trial(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory);

  /// Perform update action.
  virtual void update(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory);

  /// Perform write action.
  virtual std::string write(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory);

  // Access to factory of Analyze objects.
  virtual const std::vector<std::shared_ptr<Analyze> >& analyzers() const;
  virtual const Analyze * analyze(const int index) const;

  // serialization
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Analyze> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Analyze> >& deserialize_map();
  std::shared_ptr<Analyze> deserialize(std::istream& istr);
  explicit Analyze(std::istream& istr) : Stepper(istr) {}
  virtual ~Analyze() {}
  std::string class_name() const override { return std::string("Analyze"); }
};

/**
  This Analyze does not perform updates.
 */
class AnalyzeWriteOnly : public Analyze {
 public:
  /**
    args:
    - steps_per: write every this many steps
   */
  explicit AnalyzeWriteOnly(const argtype &args = argtype());

  void set_steps_per_update(const int steps) override;

  void set_steps_per(const int steps) { set_steps_per_write(steps); }

  explicit AnalyzeWriteOnly(std::istream& istr) : Analyze(istr) {}
};

/**
  This Analyze does not perform writes.
 */
class AnalyzeUpdateOnly : public Analyze {
 public:
  /**
    args:
    - steps_per: update every this many steps
   */
  explicit AnalyzeUpdateOnly(const argtype &args = argtype());

  void set_steps_per_write(const int steps) override;

  void set_steps_per(const int steps) { set_steps_per_update(steps); }

  explicit AnalyzeUpdateOnly(std::istream& istr) : Analyze(istr) {}
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ANALYZE_H_
