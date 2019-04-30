
#ifndef FEASST_CORE_ANALYZE_FACTORY_H_
#define FEASST_CORE_ANALYZE_FACTORY_H_

#include "core/include/analyze.h"

namespace feasst {

class AnalyzeFactory : public Analyze {
 public:
  AnalyzeFactory() : Analyze() {}

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

  // std::shared_ptr<Analyze> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  AnalyzeFactory(std::istream& istr);

 private:
  std::string class_name_ = "AnalyzeFactory";
  std::vector<std::shared_ptr<Analyze> > analyzers_;
};

}  // namespace feasst

#endif  // FEASST_CORE_ANALYZE_FACTORY_H_
