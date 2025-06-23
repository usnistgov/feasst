#ifndef FEASST_MONTE_CARLO_TRIAL_ADD_REMOVE_H_
#define FEASST_MONTE_CARLO_TRIAL_ADD_REMOVE_H_

#include <map>
#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
 * Half of the time, attempt to add a particle of given type.
 * The other half of the time, attempt to remove a particle of that type.
 */
class TrialAddRemove : public Trial {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - TrialSelectParticle arguments.
    - TrialStage arguments.
   */
  explicit TrialAddRemove(argtype args = argtype());
  explicit TrialAddRemove(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAddRemove>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAddRemove>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddRemove(std::istream& istr);
  virtual ~TrialAddRemove();
  //@}
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_ADD_REMOVE_H_
