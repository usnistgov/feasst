
#ifndef FEASST_CORE_ANALYZE_FACTORY_H_
#define FEASST_CORE_ANALYZE_FACTORY_H_

#include "core/include/analyze.h"

namespace feasst {

class AnalyzeFactory : public Analyze {
 public:
  AnalyzeFactory(const argtype &args = argtype()) : Analyze(args) {}

  void initialize(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) override;

  void add(std::shared_ptr<Analyze> analyze) {
    analyzers_.push_back(analyze); }

  const std::vector<std::shared_ptr<Analyze> >& analyzers() const {
    return analyzers_; }

  void trial(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) override;

  std::string class_name() const override { return std::string("AnalyzeFactory"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<AnalyzeFactory>(istr); }
  void serialize(std::ostream& ostr) const override;
  AnalyzeFactory(std::istream& istr);
  virtual ~AnalyzeFactory() {}

 private:
  std::vector<std::shared_ptr<Analyze> > analyzers_;
};

inline std::shared_ptr<AnalyzeFactory> MakeAnalyzeFactory(const argtype &args = argtype()) {
  return std::make_shared<AnalyzeFactory>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_ANALYZE_FACTORY_H_
