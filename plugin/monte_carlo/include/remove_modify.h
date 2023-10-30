
#ifndef FEASST_MONTE_CARLO_REMOVE_MODIFY_H_
#define FEASST_MONTE_CARLO_REMOVE_MODIFY_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Remove a Modify.
 */
class RemoveModify : public Action {
 public:
  //@{
  /** @name Arguments
    - index: index of modify to remove, in the order added.
      If -1, do nothing. (default: -1).
    - name: remove first modify with this class name, if not empty.
      (default: empty).
    - all: remove all modifiers (default: false)
   */
  explicit RemoveModify(argtype args = argtype());
  explicit RemoveModify(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<RemoveModify>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<RemoveModify>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RemoveModify(std::istream& istr);
  virtual ~RemoveModify() {}

  //@}
 private:
  int index_;
  std::string name_;
  bool all_;
};

inline std::shared_ptr<RemoveModify> MakeRemoveModify(argtype args = argtype()) {
  return std::make_shared<RemoveModify>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_REMOVE_MODIFY_H_
