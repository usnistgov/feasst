
#ifndef FEASST_MONTE_CARLO_ANALYZE_FACTORY_H_
#define FEASST_MONTE_CARLO_ANALYZE_FACTORY_H_

#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Contains multiple Analyze objects.
 */
class AnalyzeFactory : public Analyze {
 public:
  explicit AnalyzeFactory(const argtype &args = argtype()) : Analyze(args) {}

  void initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  /// Add an Analyze object.
  void add(std::shared_ptr<Analyze> analyze) { analyzers_.push_back(analyze); }

  /// Return the number.
  int num() const { return static_cast<int>(analyzers_.size()); }

  /// Return the Analyze objects.
  const std::vector<std::shared_ptr<Analyze> >& analyzers() const override {
    return analyzers_; }

  /// Return an Analyze object by index.
  const Analyze& analyze(const int index) const override {
    return const_cast<Analyze&>(*analyzers_[index]); }

  void trial(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override {
    return std::string("AnalyzeFactory"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<AnalyzeFactory>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit AnalyzeFactory(std::istream& istr);
  virtual ~AnalyzeFactory() {}

 private:
  std::vector<std::shared_ptr<Analyze> > analyzers_;

  void trial_(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory,
    const int index);
};

inline std::shared_ptr<AnalyzeFactory> MakeAnalyzeFactory(
    const argtype &args = argtype()) {
  return std::make_shared<AnalyzeFactory>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ANALYZE_FACTORY_H_
