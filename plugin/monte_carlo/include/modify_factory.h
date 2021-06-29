
#ifndef FEASST_MONTE_CARLO_MODIFY_FACTORY_H_
#define FEASST_MONTE_CARLO_MODIFY_FACTORY_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Contains multiple Modify objects.
 */
class ModifyFactory : public Modify {
 public:
  explicit ModifyFactory(argtype args = argtype()) : Modify(&args) {}

  /// Add a Modify object.
  void add(std::shared_ptr<Modify> modify) { modifiers_.push_back(modify); }

  /// Remove a Modify
  void remove(const int index) { modifiers_.erase(modifiers_.begin() + index); }

  /// Return the number.
  int num() const { return static_cast<int>(modifiers_.size()); }

  /// Return the Modify objects.
  const std::vector<std::shared_ptr<Modify> >& modifiers() const override {
    return modifiers_; }

  /// Return a Modify object by index.
  const Modify& modify(const int index) const override {
    return const_cast<Modify&>(*modifiers_[index]); }

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

  void trial_(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory,
    const int index);
};

inline std::shared_ptr<ModifyFactory> MakeModifyFactory(
    argtype args = argtype()) {
  return std::make_shared<ModifyFactory>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_MODIFY_FACTORY_H_
