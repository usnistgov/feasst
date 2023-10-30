
#ifndef FEASST_STEPPERS_CRITERIA_WRITER_H_
#define FEASST_STEPPERS_CRITERIA_WRITER_H_

#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Periodically write Criteria.
 */
class CriteriaWriter : public AnalyzeWriteOnly {
 public:
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit CriteriaWriter(argtype args = argtype());
  explicit CriteriaWriter(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override {
    return std::string("CriteriaWriter"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<CriteriaWriter>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<CriteriaWriter>(args); }
  CriteriaWriter(std::istream& istr);
  //@}
};

inline std::shared_ptr<CriteriaWriter> MakeCriteriaWriter(
    argtype args = argtype()) {
  return std::make_shared<CriteriaWriter>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CRITERIA_WRITER_H_
