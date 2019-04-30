
#ifndef FEASST_CORE_MODIFY_FACTORY_H_
#define FEASST_CORE_MODIFY_FACTORY_H_

#include "core/include/modify.h"

namespace feasst {

class ModifyFactory : public Modify {
 public:
  ModifyFactory() {}

  void initialize(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) override {
    for (const std::shared_ptr<Modify> modify : modifiers_) {
      modify->initialize(criteria, system, trial_factory);
    }
  }

  void add(std::shared_ptr<Modify> modify) {
    modifiers_.push_back(modify);
  }

  const std::vector<std::shared_ptr<Modify> >& modifiers() const {
    return modifiers_; }

  void trial(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) override {
    for (const std::shared_ptr<Modify> modify : modifiers_) {
      modify->trial(criteria, system, trial_factory);
    }
  }

  void serialize(std::ostream& ostr) const override;
  ModifyFactory(std::istream& istr);

 private:
  std::string class_name_ = "ModifyFactory";
  std::vector<std::shared_ptr<Modify> > modifiers_;
};

}  // namespace feasst

#endif  // FEASST_CORE_MODIFY_FACTORY_H_
