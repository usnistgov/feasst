
#ifndef FEASST_ACTIONS_WAIT_H_
#define FEASST_ACTIONS_WAIT_H_

#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include "monte_carlo/include/action.h"

namespace feasst {

class Analyze;
class Modify;
class Trial;
class TrialFactoryNamed;

/**
  Perform a number of trials.
 */
class Wait : public Action {
 public:
  //@{
  /** @name Arguments
    - until_file_exists: wait until the given file name exists.
      If multiple files are required, then use the following syntax:
      "prefix[number]suffix" where number is the number of files and [number]
      is replaced with the interger indices of [0, number-1] in increments of
      one.
    - seconds_per_check: wait this many seconds before checking (default: 5).
   */
  explicit Wait(argtype args = argtype());
  explicit Wait(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  const std::string& prefix() const { return prefix_; }
  const std::string& suffix() const { return suffix_; }
  int number() const { return number_; }

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<Wait>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<Wait>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Wait(std::istream& istr);
  virtual ~Wait() {}

  //@}
 private:
  std::string until_file_exists_, prefix_, suffix_;
  double seconds_per_check_;
  int number_ = 0;
};

}  // namespace feasst

#endif  // FEASST_ACTIONS_WAIT_H_
