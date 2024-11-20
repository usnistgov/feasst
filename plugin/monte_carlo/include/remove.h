
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
    - name[i]: remove first class with this class name, if not empty.
      (default: empty).
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one remove, the "[i]" is optional.
    - name_contains[i]: same "name[i]" above, except the entire name does not
      have to match exactly and all matches are removed (not just the first).
      If any part of the given characters match then remove, if not empty.
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
