
#ifndef FEASST_MONTE_CARLO_REMOVE_H_
#define FEASST_MONTE_CARLO_REMOVE_H_

#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Remove an Analyze, Modify or Trial based on the name.
 */
class Remove : public Action {
 public:
  //@{
  /** @name Arguments
    - name: remove first class with this class name, if not empty.
      (default: empty).
      Multiple names can be provided as comma-separated values.
    - name_contains: same as name above, except the entire name does not
      have to match exactly and all matches are removed (not just the first).
      Multiple name_contains can be provided as comma-separated values.
    - all_trials: if true (default: false), remove all trials.
    - all_analyzers: if true (default: false), remove all analyzers.
    - all_modifiers: if true (default: false), remove all modifiers.
   */
  explicit Remove(argtype args = argtype());
  explicit Remove(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<Remove>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<Remove>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Remove(std::istream& istr);
  virtual ~Remove() {}

  //@}
 private:
  std::vector<std::string> names_, name_contains_;
  bool all_trials_, all_analyzers_, all_modifiers_;
};

inline std::shared_ptr<Remove> MakeRemove(argtype args = argtype()) {
  return std::make_shared<Remove>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_REMOVE_H_
