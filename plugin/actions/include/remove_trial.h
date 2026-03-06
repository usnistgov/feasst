
#ifndef FEASST_MONTE_CARLO_REMOVE_TRIAL_H_
#define FEASST_MONTE_CARLO_REMOVE_TRIAL_H_

#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Remove a Trial.
 */
class RemoveTrial : public Action {
 public:
  //@{
  /** @name Arguments
    - index: index of trial to remove, in order of trials added.
      If -1, do nothing. (default: -1).
    - name: remove first trial with this class name, if not empty.
      (default: empty).
    - all: remove all trials (default: false)
    - name_contains: remove all trials with this in the class name, if not
      empty (default: empty).
      Similar to the name that appears in Log, if the class_name is Trial,
      use the description instead (e.g., TrialGrowadd in TrialGrow).
   */
  explicit RemoveTrial(argtype args = argtype());
  explicit RemoveTrial(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<RemoveTrial>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<RemoveTrial>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RemoveTrial(std::istream& istr);
  virtual ~RemoveTrial() {}

  //@}
 private:
  int index_;
  std::string name_;
  std::string name_contains_;
  bool all_;
};

inline std::shared_ptr<RemoveTrial> MakeRemoveTrial(argtype args = argtype()) {
  return std::make_shared<RemoveTrial>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_REMOVE_TRIAL_H_
