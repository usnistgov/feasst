
#ifndef FEASST_STEPPERS_CRITERIA_UPDATER_H_
#define FEASST_STEPPERS_CRITERIA_UPDATER_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Periodically write Criteria.
 */
class CriteriaUpdater : public ModifyUpdateOnly {
 public:
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit CriteriaUpdater(argtype args = argtype());
  explicit CriteriaUpdater(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) override { criteria->update(); }

  // serialize
  std::string class_name() const override {
    return std::string("CriteriaUpdater"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<CriteriaUpdater>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<CriteriaUpdater>(args); }
  CriteriaUpdater(std::istream& istr);
  //@}
};

inline std::shared_ptr<CriteriaUpdater> MakeCriteriaUpdater(
    argtype args = argtype()) {
  return std::make_shared<CriteriaUpdater>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CRITERIA_UPDATER_H_
