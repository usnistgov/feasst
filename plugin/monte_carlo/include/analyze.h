
#ifndef FEASST_MONTE_CARLO_ANALYZE_H_
#define FEASST_MONTE_CARLO_ANALYZE_H_

#include <vector>
#include <memory>
#include <string>
#include <map>
#include "monte_carlo/include/stepper.h"

namespace feasst {

class Configuration;
class MonteCarlo;
class TrialFactory;

/**
  Perform an Action every so many attempted trials that will not change the
  System, Criteria or Trials.
 */
class Analyze : public Stepper {
 public:
  Analyze() : Stepper() {}
  explicit Analyze(argtype * args) : Stepper(args) {}

  /// Initialize and precompute before trials.
  void initialize(MonteCarlo * mc) override { Stepper::initialize(mc); }

  /// Check every trial if action is to be performed.
  virtual void trial(const MonteCarlo& mc);

  /// Perform update action.
  virtual void update(const MonteCarlo& mc);

  /// Perform write action.
  virtual std::string write(const MonteCarlo& mc);
  virtual void write_to_file(const MonteCarlo& mc);

  // Access to factory of Analyze objects.
  virtual const std::vector<std::shared_ptr<Analyze> >& analyzers() const;
  virtual const Analyze& analyze(const int index) const;
  virtual Analyze * get_analyze(const int index);

  // serialization
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Analyze> create(std::istream& istr) const;
  virtual std::shared_ptr<Analyze> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<Analyze> >& deserialize_map();
  std::shared_ptr<Analyze> deserialize(std::istream& istr);
  std::shared_ptr<Analyze> factory(const std::string name, argtype * args);
  explicit Analyze(std::istream& istr) : Stepper(istr) {}
  virtual ~Analyze() {}
  std::string class_name() const override { return std::string("Analyze"); }

  // HWH only used by AnalyzeFactory
  void check_update_(const MonteCarlo& mc);
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ANALYZE_H_
