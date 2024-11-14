
#ifndef FEASST_GIBBS_COPY_PREVIOUS_LINE_H_
#define FEASST_GIBBS_COPY_PREVIOUS_LINE_H_

#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Copy the next line that was parsed by MonteCarlo.
 */
class CopyNextLine : public Action {
 public:
  //@{
  /** @name Arguments
    - replace[i]: replace this argument with the following with[i].
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one replace, the "[i]" is optional.
    - with[i]: this is the value to replace the argument with.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one with, the "[i]" is optional.
   */
  explicit CopyNextLine(argtype args = argtype());
  explicit CopyNextLine(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<CopyNextLine>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<CopyNextLine>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit CopyNextLine(std::istream& istr);
  virtual ~CopyNextLine() {}

  //@}
 private:
  std::vector<std::vector<std::string> > replace_;
};

}  // namespace feasst

#endif  // FEASST_GIBBS_COPY_PREVIOUS_LINE_H_
