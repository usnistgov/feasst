
#ifndef FEASST_MONTE_CARLO_REMOVE_ANALYZE_H_
#define FEASST_MONTE_CARLO_REMOVE_ANALYZE_H_

#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Remove an Analyze.
 */
class RemoveAnalyze : public Action {
 public:
  //@{
  /** @name Arguments
    - index: index of analyze to remove, in the order added.
      If -1, do nothing. (default: -1).
    - name: remove first analyze with this class name, if not empty.
      (default: empty).
    - all: remove all analyzers (default: false)
   */
  explicit RemoveAnalyze(argtype args = argtype());
  explicit RemoveAnalyze(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<RemoveAnalyze>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<RemoveAnalyze>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RemoveAnalyze(std::istream& istr);
  virtual ~RemoveAnalyze() {}

  //@}
 private:
  int index_;
  std::string name_;
  bool all_;
};

inline std::shared_ptr<RemoveAnalyze> MakeRemoveAnalyze(
    argtype args = argtype()) {
  return std::make_shared<RemoveAnalyze>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_REMOVE_ANALYZE_H_
