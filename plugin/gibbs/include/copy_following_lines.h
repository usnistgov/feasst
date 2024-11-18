
#ifndef FEASST_GIBBS_COPY_FOLLOWING_LINES_H_
#define FEASST_GIBBS_COPY_FOLLOWING_LINES_H_

#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Copy the following lines for each configuration_index until EndCopy is reached.
 */
class CopyFollowingLines : public Action {
 public:
  //@{
  /** @name Arguments
    - for_num_configurations: number of configurations to copy (default: 1).
    - replace[i]: replace this argument with the following with[i].
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one replace, the "[i]" is optional.
    - with[i]: this is the value to replace the argument with.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one with, the "[i]" is optional.
   */
  explicit CopyFollowingLines(argtype args = argtype());
  explicit CopyFollowingLines(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<CopyFollowingLines>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<CopyFollowingLines>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit CopyFollowingLines(std::istream& istr);
  virtual ~CopyFollowingLines() {}

  //@}
 private:
  std::vector<std::vector<std::string> > replace_;
  int for_num_configurations_;
};

}  // namespace feasst

#endif  // FEASST_GIBBS_COPY_FOLLOWING_LINES_H_
