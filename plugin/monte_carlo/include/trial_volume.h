#ifndef FEASST_MONTE_CARLO_TRIAL_VOLUME_H_
#define FEASST_MONTE_CARLO_TRIAL_VOLUME_H_

#include <map>
#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Attempt to change the volume.
  Not currently fully implemented or tested for many cases, including Ewald
  or VisitModelCell, etc.
 */
class TrialVolume : public Trial {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - PerturbVolume arguments.
    - Tunable arguments.
   */
  explicit TrialVolume(argtype args = argtype());
  explicit TrialVolume(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialVolume>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialVolume>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialVolume(std::istream& istr);
  virtual ~TrialVolume() {}
  //@}
};

inline std::shared_ptr<TrialVolume> MakeTrialVolume(argtype args = argtype()) {
  return std::make_shared<TrialVolume>(args); }

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_VOLUME_H_
