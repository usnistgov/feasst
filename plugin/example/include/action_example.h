
#ifndef FEASST_EXAMPLE_ACTION_EXAMPLE_H_
#define FEASST_EXAMPLE_ACTION_EXAMPLE_H_

#include <string>
#include <memory>
#include <map>
#include "monte_carlo/include/action.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Add an action to FEASST by using this file as a template and instruction set.
  Follow the same steps detailed in /feasst/plugin/example/README.rst.
  In summary, copy action_example.[h/cpp] to new_name.[h/cpp], replace
  ActionExample with NewName, then replace ACTION_EXAMPLE with NEW_NAME.

  An Action is serializable, meaning that a restarted simulation will still
  complete the action in the order specified in the text file.
  Unlike Stepper/Analyze/Modify, Action happens only once per instance.

  In this example, an Analyze/Modify of given name is forced to write to file.
  This is convenient for example printing the very last update to FlatHistogram.
  This is a duplicate of WriteStepper class.
 */
class ActionExample : public Action {
 public:
  //@{
  /** @name Arguments
    - analyze_name: The Analyze::class_name() to print to file.
      If empty, do nothing (default: empty).
    - modify_name: The Modify::class_name() to print to file (default: empty).
      If empty, do nothing (default: empty).
   */
  explicit ActionExample(argtype args = argtype());
  explicit ActionExample(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<ActionExample>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<ActionExample>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ActionExample(std::istream& istr);
  virtual ~ActionExample() {}

  //@}
 private:
  std::string analyze_name_;
  std::string modify_name_;
};

}  // namespace feasst

#endif  // FEASST_EXAMPLE_ACTION_EXAMPLE_H_
