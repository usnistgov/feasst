
#ifndef FEASST_MONTE_CARLO_TRIAL_MOVE_H_
#define FEASST_MONTE_CARLO_TRIAL_MOVE_H_

#include <map>
#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

class PerturbMove;
class TrialSelect;

typedef std::map<std::string, std::string> argtype;

/// Attempt to rigidly move a selection in a Trial in one stage.
std::shared_ptr<Trial> MakeTrialMove(
    std::shared_ptr<TrialSelect> select,
    std::shared_ptr<PerturbMove> perturb,
    const std::string& description,
    argtype * args);

class TrialMove : public Trial {
 public:
  TrialMove(std::shared_ptr<TrialSelect> select,
            std::shared_ptr<PerturbMove> perturb,
            argtype * args);
  explicit TrialMove(std::istream& istr);
  virtual ~TrialMove() {}
 protected:
  void serialize_trial_move_(std::ostream& ostr) const;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_MOVE_H_
