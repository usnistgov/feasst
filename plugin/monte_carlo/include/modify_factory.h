
#ifndef FEASST_MONTE_CARLO_MODIFY_FACTORY_H_
#define FEASST_MONTE_CARLO_MODIFY_FACTORY_H_

#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/modify.h"

namespace feasst {

class TimerRDTSC;

/**
  Contains multiple Modify objects.
 */
class ModifyFactory : public Modify {
 public:
  explicit ModifyFactory(argtype args = argtype());

  /// Add a Modify object.
  void add(std::shared_ptr<Modify> modify);

  /// Remove a Modify
  void remove(const int index);

  /// Return the number.
  int num() const { return static_cast<int>(modifiers_.size()); }

  /// Return the Modify objects.
  const std::vector<std::shared_ptr<Modify> >& modifiers() const override {
    return modifiers_; }

  /// Return a Modify object by index.
  const Modify& modify(const int index) const override {
    return const_cast<Modify&>(*modifiers_[index]); }

  void initialize(MonteCarlo * mc) override;

  /// Write all Modify immediately.
  void write_to_file(MonteCarlo * mc) override;

  /// For use with CollectionMatrixSplice, transfer multistate between threads.
  void adjust_bounds(const bool adjusted_up, const std::vector<int>& states,
    ModifyFactory * analyze_factory);

  void trial(MonteCarlo * mc) override;

  Modify * get_modify(const int index) override {
    return modifiers_[index].get(); }

  /// Set the timer
  void set_timer();

  /// Return timer
  const TimerRDTSC * const timer() const { return timer_.get(); }

  void synchronize_(const Modify& modify) override;

  std::string class_name() const override {
    return std::string("ModifyFactory"); }
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<ModifyFactory>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit ModifyFactory(std::istream& istr);
  virtual ~ModifyFactory();

 private:
  std::vector<std::shared_ptr<Modify> > modifiers_;
  std::unique_ptr<TimerRDTSC> timer_;

  void trial_(MonteCarlo * mc, const int index);
};

inline std::shared_ptr<ModifyFactory> MakeModifyFactory(
    argtype args = argtype()) {
  return std::make_shared<ModifyFactory>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_MODIFY_FACTORY_H_
