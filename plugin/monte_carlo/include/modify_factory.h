
#ifndef FEASST_MONTE_CARLO_MODIFY_FACTORY_H_
#define FEASST_MONTE_CARLO_MODIFY_FACTORY_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Contains multiple Modify objects.
 */
class ModifyFactory : public Modify {
 public:
  explicit ModifyFactory(const argtype &args = argtype()) : Modify(args) {}

  /// Add a Modify object.
  void add(std::shared_ptr<Modify> modify) { modifiers_.push_back(modify); }

  /// Return the Modify objects.
  const std::vector<std::shared_ptr<Modify> >& modifiers() const {
    return modifiers_; }

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  void trial(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  std::string class_name() const override {
    return std::string("ModifyFactory"); }
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<ModifyFactory>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit ModifyFactory(std::istream& istr);
  virtual ~ModifyFactory() {}

 private:
  std::vector<std::shared_ptr<Modify> > modifiers_;
};

inline std::shared_ptr<ModifyFactory> MakeModifyFactory(
    const argtype &args = argtype()) {
  return std::make_shared<ModifyFactory>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_MODIFY_FACTORY_H_
