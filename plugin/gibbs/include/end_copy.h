
#ifndef FEASST_GIBBS_END_COPY_H_
#define FEASST_GIBBS_END_COPY_H_

#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Copy the following lines for each configuration_index until EndCopy is reached.
 */
class EndCopy : public Action {
 public:
  //@{
  /** @name Arguments
   */
  explicit EndCopy(argtype args = argtype());
  explicit EndCopy(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<EndCopy>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<EndCopy>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit EndCopy(std::istream& istr);
  virtual ~EndCopy() {}

  //@}
 private:
};

}  // namespace feasst

#endif  // FEASST_GIBBS_END_COPY_H_
